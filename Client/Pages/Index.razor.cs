﻿using Microsoft.JSInterop;
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
            await JSRuntime.InvokeVoidAsync("QuantumInterop.initializeThreeJs", 
                new { x = BoxX, y = BoxY, z = BoxZ }, SpaceStep);
            await RestartSimulation();
        }
    }

    private CancellationTokenSource _cancellationTokenSource = new();
    private async Task RestartSimulation()
    {
        _cancellationTokenSource.Cancel();
        _cancellationTokenSource.Dispose();
        _cancellationTokenSource = new CancellationTokenSource();

        var newSystem = await ResetSimAndGraphics();
        await StartSimulation(newSystem, _cancellationTokenSource.Token);
    }

    private async Task<QuantumSystem> ResetSimAndGraphics()
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

        await Update3DProbabilityDisplay(newSystem);
        StateHasChanged();
        return newSystem;
    }

    private async Task StartSimulation(QuantumSystem quantumSystem, CancellationToken token)
    {
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
                        Update3DProbabilityDisplay(quantumSystem), "UpdateProbabilityData(quantumSystem)");
                }, token);

                if (token.IsCancellationRequested) break;
                
                await Profiling.RunWithClockingLogAsync(() => 
                    DelayUntilNextFrame(token), "DelayUntilNextFrame(token)");

                await Profiling.RunWithClockingLogAsync(() => 
                    InvokeAsync(StateHasChanged), "InvokeAsync(StateHasChanged)");
            }
        }
    }

    private DateTime _lastFrameTime = DateTime.UtcNow;
    private async Task DelayUntilNextFrame(CancellationToken token)
    {
        var currentFrameTime = DateTime.UtcNow;
        var elapsed = currentFrameTime - _lastFrameTime;

        if (elapsed.TotalSeconds > 0)
            FrameRate = 1.0f / (float)elapsed.TotalSeconds;
        _lastFrameTime = currentFrameTime;
        FrameCount++;

        var difference = (1000 / MaxFramesPerSecond) - (int)elapsed.TotalMilliseconds;
        if (difference > 0)
            await Task.Delay(difference, token);
    }

    private async Task Update3DProbabilityDisplay(QuantumSystem quantumSystem)
    {
        var filteredData = quantumSystem.GetSignificantlyChangedProbability();
        if (filteredData.Count > 0)
            await JSRuntime.InvokeVoidAsync("QuantumInterop.updateThreeJsScene", filteredData);
    }

    private void UpdateEnergy(QuantumSystem quantumSystem)
    {
        CurrentTotalEnergy = quantumSystem.CalculateTotalEnergy();
        if (OriginalTotalEnergy == 0) OriginalTotalEnergy = CurrentTotalEnergy;
    }

    private void TogglePause() => Paused = !Paused;
    private void ToggleManualInputs() => ShowManualInputs = !ShowManualInputs;
    private void ToggleControlsVisibility() => AreControlsVisible = !AreControlsVisible;
    private string GetControlPanelClass() => AreControlsVisible ? "expanded" : "collapsed";

    private async Task InstallApp()
    {
        var success = await JSRuntime.InvokeAsync<bool>("showPWAInstallPrompt");
        Console.WriteLine(
            $"Install prompt was {(success ? "" : "not")} accepted." +
            $"{(success ? "" : " Event to trigger it might not have been caught.")}");
    }
}