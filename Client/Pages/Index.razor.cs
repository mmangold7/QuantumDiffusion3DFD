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

    private bool _isSimulationRunning;
    private bool _areControlsVisible = true;
    private bool _showManualInputs;
    private int _frameCount;
    private float _frameRate;
    private DateTime _lastFrameTime = DateTime.UtcNow;
    private float[]? _previousProbabilityData;
    private CancellationTokenSource _cancellationTokenSource = new();

    private bool Paused { get; set; }
    private int MaxFramesPerSecond { get; set; } = 30;
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
            await JSRuntime.InvokeVoidAsync("QuantumInterop.initializeThreeJs", 
                new { x = BoxX, y = BoxY, z = BoxZ }, SpaceStep);
            await RestartSimulation();
        }
    }

    private async Task InstallApp()
    {
        var success = await JSRuntime.InvokeAsync<bool>("showPWAInstallPrompt");
        Console.WriteLine(
            $"Install prompt was {(success ? "" : "not")} accepted." +
            $"{(success ? "" : " Event to trigger it might not have been caught.")}");
    }

    private async Task<QuantumSystem> ResetSimAndGraphics()
    {
        _frameCount = 0;
        _frameRate = 0;
        OriginalTotalEnergy = 0;
        CurrentTotalEnergy = 0;
        _previousProbabilityData = null;

        var newSystem = new QuantumSystem(
            BoxX, BoxY, BoxZ, BoundaryType, TimeStep, SpaceStep, Mass, Hbar);
        newSystem.InitializeGaussianPacket(
            ParticleXPosition, ParticleYPosition, ParticleZPosition, 
            RadiusSigma, 
            XMomentumWavenumber, YMomentumWavenumber, ZMomentumWavenumber);

        await UpdateProbabilityData(newSystem);
        StateHasChanged();
        return newSystem;
    }

    private async Task StartSimulation(QuantumSystem quantumSystem, CancellationToken token)
    {
        _isSimulationRunning = true;

        while (!token.IsCancellationRequested)
        {
            if (Paused) await Task.Delay(100, token); // Small delay to reduce CPU usage
            else
            {
                await Task.Run(async () =>
                {
                    Profiling.RunWithClockingLog(quantumSystem.ApplySingleTimeEvolutionStepEuler);

                    Profiling.RunWithClockingLog(() => 
                        UpdateEnergy(quantumSystem), "UpdateEnergy(quantumSystem)"); 

                    await Profiling.RunWithClockingLogAsync(() => 
                        UpdateProbabilityData(quantumSystem), "UpdateProbabilityData(quantumSystem)");
                }, token);

                if (token.IsCancellationRequested) break;
                
                await Profiling.RunWithClockingLogAsync(() => 
                    DelayUntilNextFrame(token), "DelayUntilNextFrame(token)");

                await Profiling.RunWithClockingLogAsync(() => 
                    InvokeAsync(StateHasChanged), "InvokeAsync(StateHasChanged)");
            }
        }

        _isSimulationRunning = false;
    }

    private async Task UpdateProbabilityData(QuantumSystem quantumSystem)
    {
        var newProbabilityData = quantumSystem.GetNormalizedProbabilityData();
        var filteredData = FilterSignificantUpdates(newProbabilityData);
        await Update3DDisplay(filteredData);
    }

    private void UpdateEnergy(QuantumSystem quantumSystem)
    {
        CurrentTotalEnergy = quantumSystem.CalculateTotalEnergy();
        if (OriginalTotalEnergy == 0) OriginalTotalEnergy = CurrentTotalEnergy;
    }

    private async Task DelayUntilNextFrame(CancellationToken token)
    {
        var currentFrameTime = DateTime.UtcNow;
        var elapsed = currentFrameTime - _lastFrameTime;

        if (elapsed.TotalSeconds > 0)
            _frameRate = 1.0f / (float)elapsed.TotalSeconds;
        _lastFrameTime = currentFrameTime;
        _frameCount++;

        var difference = (1000 / MaxFramesPerSecond) - (int)elapsed.TotalMilliseconds;
        if (difference > 0)
            await Task.Delay(difference, token);
    }

    public async Task Update3DDisplay(List<object> probabilityData)
    {
        if (probabilityData.Count > 0)
            await JSRuntime.InvokeVoidAsync("QuantumInterop.updateThreeJsScene", probabilityData);
    }

    // Prevents lag from "updating" all of the cubes in three.js that have barely changed in value
    // Ideally threshold would be calculated based on whether update would be perceivable/detectable by user (because of significantly different opacity and/or color)
    private List<object> FilterSignificantUpdates(float[] probabilityData)
    {
        var updatedData = new List<object>();
        float maxProbability = 1;
        var updateThreshold = maxProbability * 0.1f;
        var opacityScale = 0.75f;

        for (int i = 0; i < probabilityData.Length; i++)
        {
            var newProbability = probabilityData[i];
            if (_previousProbabilityData == null ||
                Math.Abs(newProbability - _previousProbabilityData[i]) > updateThreshold)
            {
                var color = GraphicsExtensions.InterpolateColor(newProbability);
                var opacity = GraphicsExtensions.SigmoidOpacity(newProbability * opacityScale);
                updatedData.Add(new { index = i, color, opacity });
                if (_previousProbabilityData != null) _previousProbabilityData[i] = newProbability;
            }
        }

        _previousProbabilityData ??= probabilityData;
        return updatedData;
    }

    private async Task RestartSimulation()
    {
        if (_isSimulationRunning)
        {
            _cancellationTokenSource.Cancel();
            _cancellationTokenSource.Dispose();
            _cancellationTokenSource = new CancellationTokenSource();
        }

        var newSystem = await ResetSimAndGraphics();
        await StartSimulation(newSystem, _cancellationTokenSource.Token);
    }

    private void TogglePause() => Paused = !Paused;
    private void ToggleManualInputs() => _showManualInputs = !_showManualInputs;
    private void ToggleControlsVisibility() => _areControlsVisible = !_areControlsVisible;
    private string GetControlPanelClass() => _areControlsVisible ? "expanded" : "collapsed";
}