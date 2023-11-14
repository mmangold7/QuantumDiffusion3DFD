using QuantumDiffusion3DFD.Shared;
using System.Numerics;

namespace QuantumDiffusion3DFD.Client.Pages;

public partial class Index
{
    private QuantumSystem quantumSystem;
    private double simulationTime;
    private Complex wavefunctionValue;
    private int frameCount = 0;
    private double frameRate = 0.0;

    public bool Paused { get; set; }
    public bool RunningSimulation { get; set; }
    public int SimulationDelayMilliseconds { get; set; } = 33;

    protected override async Task OnInitializedAsync()
    {
        quantumSystem = new QuantumSystem((10, 10, 10), 0.001, 0.1, BoundaryType.Reflective);
        quantumSystem.InitializeWavefunction(new Complex(1.0, 0.0));

        RunningSimulation = true;
        await QuantumSimulationLoop(default);
    }

    private void TogglePause() => Paused = !Paused;

    private async Task QuantumSimulationLoop(CancellationToken token)
    {
        if (quantumSystem == null) return;
        
        DateTime lastFrameTime = DateTime.UtcNow;
        while (RunningSimulation)
        {
            DateTime currentFrameTime = DateTime.UtcNow;
            double elapsedSeconds = (currentFrameTime - lastFrameTime).TotalSeconds;
            frameRate = 1.0 / elapsedSeconds;
            lastFrameTime = currentFrameTime;

            frameCount++;

            //ProcessUserInputs();

            if (!Paused)
            {
                quantumSystem.ApplySingleTimeEvolutionStep();
                UpdateDisplay();
            }

            await Task.Delay(SimulationDelayMilliseconds, token);
        }
    }

    private void UpdateDisplay()
    {
        int x = 5;
        int y = 5;
        int z = 5;
        wavefunctionValue = quantumSystem.GetWavefunctionValue(x, y, z);
        StateHasChanged();
    }
}