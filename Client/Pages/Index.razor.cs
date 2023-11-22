using Microsoft.JSInterop;
using QuantumDiffusion3DFD.Shared;

namespace QuantumDiffusion3DFD.Client.Pages;

public partial class Index
{
    private const int CubeDimensionSize = 10;
    private const int BoxX = CubeDimensionSize;
    private const int BoxY = CubeDimensionSize;
    private const int BoxZ = CubeDimensionSize;
    private const float SpaceStep = 0.1f;

    private bool AreControlsVisible { get; set; } = true;
    private bool ShowManualInputs { get; set; }
    private bool Paused { get; set; }
    private int MaxFramesPerSecond { get; set; } = 30;
    private int FrameCount { get; set; }
    private float FrameRate { get; set; }
    private float TimeStep { get; set; } = 0.01f;
    private float ParticleXPosition { get; set; } = BoxX / 2.0f;
    private float ParticleYPosition { get; set; } = BoxY / 2.0f;
    private float ParticleZPosition { get; set; } = BoxZ / 2.0f;
    private float RadiusSigma { get; set; } = 2;
    private float XMomentumWavenumber { get; set; } = 1;
    private float YMomentumWavenumber { get; set; } = 1;
    private float ZMomentumWavenumber { get; set; } = 1;
    private float Mass { get; set; } = 1000;
    private float Hbar { get; set; } = 10;
    private float OriginalTotalEnergy { get; set; }
    private float CurrentTotalEnergy { get; set; }
    private BoundaryType BoundaryType { get; set; } = BoundaryType.Reflective;

    protected override async Task OnAfterRenderAsync(bool firstRender)
    {
        if (firstRender)
        {
            await Init3dDisplay();
            await RestartSimulation();
        }
    }

    private CancellationTokenSource _simLoopCancel = new();
    private async Task RestartSimulation()
    {
        _simLoopCancel.Cancel();
        _simLoopCancel.Dispose();
        _simLoopCancel = new CancellationTokenSource();

        var newSystem = await ResetSimAndGraphics(_simLoopCancel.Token);
        await StartSimulation(newSystem, _simLoopCancel.Token);
    }

    private async Task<QuantumSystem> ResetSimAndGraphics(CancellationToken simLoopCancelToken)
    {
        FrameCount = 0;
        FrameRate = 0;
        OriginalTotalEnergy = 0;
        CurrentTotalEnergy = 0;

        var newSystem = new QuantumSystem(
            BoxX, BoxY, BoxZ, BoundaryType, TimeStep, SpaceStep, Mass, Hbar);
        newSystem.InitializeGaussianPacket(
            ParticleXPosition, ParticleYPosition, ParticleZPosition, 
            RadiusSigma, 
            XMomentumWavenumber, YMomentumWavenumber, ZMomentumWavenumber);

        await UpdateUi(newSystem, true, simLoopCancelToken);
        StateHasChanged();
        return newSystem;
    }

    private async Task StartSimulation(QuantumSystem quantumSystem, CancellationToken simLoopCancelToken)
    {
        while (!simLoopCancelToken.IsCancellationRequested)
        {
            if (Paused) await Task.Delay(100, simLoopCancelToken); // Small delay to reduce CPU usage
            else await Task.Run(async () => 
                await UpdateUi(quantumSystem, false, simLoopCancelToken), simLoopCancelToken);
            //if (simLoopCancelToken.IsCancellationRequested) break;
        }
    }

    private async Task UpdateUi(QuantumSystem newSystem, bool useAllData, CancellationToken simLoopCancelToken)
    {
        var newState = newSystem.UpdateSimulation(!useAllData);
        UpdateEnergyDisplay(newState.CurrentTotalEnergy);
        var filteredData = newState.ProbabilityData;
        if (filteredData.Count > 0) await Update3dDisplay(filteredData);
        await DelayUntilNextRequestedFrame(simLoopCancelToken);
        await InvokeAsync(StateHasChanged);
    }

    private DateTime _lastFrameTime = DateTime.UtcNow;
    private async Task DelayUntilNextRequestedFrame(CancellationToken simLoopCancelToken)
    {
        var currentFrameTime = DateTime.UtcNow;
        var elapsed = currentFrameTime - _lastFrameTime;

        if (elapsed.TotalSeconds > 0)
            FrameRate = 1.0f / (float)elapsed.TotalSeconds;
        _lastFrameTime = currentFrameTime;
        FrameCount++;

        var difference = (int)(Math.Floor(1000.0f / MaxFramesPerSecond) - (int)elapsed.TotalMilliseconds);
        if (difference > 0)
            await Task.Delay(difference, simLoopCancelToken);
    }

    private void TogglePause() => Paused = !Paused;
    private void ToggleManualInputs() => ShowManualInputs = !ShowManualInputs;
    private void ToggleControlsVisibility() => AreControlsVisible = !AreControlsVisible;
    private string GetControlPanelClass() => AreControlsVisible ? "expanded" : "collapsed";

    private void UpdateEnergyDisplay(float newEnergy)
    {
        CurrentTotalEnergy = newEnergy;
        if (OriginalTotalEnergy == 0) OriginalTotalEnergy = CurrentTotalEnergy;
    }

    private async Task Init3dDisplay() =>
        await JSRuntime.InvokeVoidAsync("QuantumInterop.initializeThreeJs",
            new { x = BoxX, y = BoxY, z = BoxZ }, SpaceStep);

    private async Task Update3dDisplay(List<object> filteredData) =>
        await JSRuntime.InvokeVoidAsync("QuantumInterop.updateThreeJsScene", filteredData);

    private async Task InstallApp()
    {
        var success = await JSRuntime.InvokeAsync<bool>("showPWAInstallPrompt");
        Console.WriteLine(
            $"Install prompt was {(success ? "" : "not")} accepted." +
            $"{(success ? "" : " Event to trigger it might not have been caught.")}");
    }
}