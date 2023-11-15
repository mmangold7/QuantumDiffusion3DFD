using QuantumDiffusion3DFD.Shared;
using System.Numerics;
using Microsoft.JSInterop;
using Microsoft.AspNetCore.Components;

namespace QuantumDiffusion3DFD.Client.Pages;

public partial class Index
{
    [Parameter]
    public float[] Data { get; set; }
    private QuantumSystem _quantumSystem;
    private Complex _wavefunctionValue;
    private int _frameCount;
    private double _frameRate;
    private (int x,int y,int z) dimensions = new(10, 10, 10);
    private double spacing = 0.1;

    public bool Paused { get; set; }
    public bool RunningSimulation { get; set; }
    public int SimulationDelayMilliseconds { get; set; } = 33;

    private void TogglePause() => Paused = !Paused;

    protected override async Task OnAfterRenderAsync(bool firstRender)
    {
        if (firstRender)
        {
            await JSRuntime.InvokeVoidAsync("initializeThreeJs", new { x = dimensions.x, y = dimensions.y, z = dimensions.z }, spacing);
        }
    }

    protected override async Task OnInitializedAsync()
    {
        SetupQuantumSystem();
        await StartQuantumSimulationLoop(default);
    }

    private void SetupQuantumSystem()
    {
        if (_quantumSystem != null) return;
        _quantumSystem = new QuantumSystem(dimensions, 0.1, spacing, BoundaryType.Reflective);
        _quantumSystem.InitializeGaussianPacket(_quantumSystem.Dimensions.x / 2, _quantumSystem.Dimensions.y / 2, _quantumSystem.Dimensions.z / 2, 1, 1, 1, 1);
    }

    private async Task StartQuantumSimulationLoop(CancellationToken token, bool paused = true)
    {
        RunningSimulation = true;
        Paused = paused;

        DateTime lastFrameTime = DateTime.UtcNow;

        while (RunningSimulation)
        {
            if (!Paused)
            {
                DateTime currentFrameTime = DateTime.UtcNow;
                double elapsedSeconds = (currentFrameTime - lastFrameTime).TotalSeconds;
                _frameRate = 1.0 / elapsedSeconds;
                lastFrameTime = currentFrameTime;

                _frameCount++;

                //ProcessUserInputs();

                _quantumSystem.ApplySingleTimeEvolutionStep();
                await UpdateDisplay();
            }

            await Task.Delay(SimulationDelayMilliseconds, token);
        }
    }

    //private void UpdateDisplay()
    //{
    //    _wavefunctionValue = _quantumSystem.GetWavefunctionValue(5, 5, 5);
    //    StateHasChanged();
    //}

    public float[] GetProbabilityData(QuantumSystem quantumSystem)
    {
        var probabilityDensity = quantumSystem.CalculateProbabilityDensity();
        var flattenedData = new List<float>();

        // Assuming the grid _quantumSystem.Dimensions in the C# simulation match the Three.js grid
        for (int x = 0; x < _quantumSystem.Dimensions.x; x++)
        {
            for (int y = 0; y < _quantumSystem.Dimensions.y; y++)
            {
                for (int z = 0; z < _quantumSystem.Dimensions.z; z++)
                {
                    // Convert complex probability to a float representation
                    // Here we just take the magnitude, but you could use other mappings
                    flattenedData.Add((float)probabilityDensity[x, y, z]);
                }
            }
        }

        // Optional: Normalize the data
        float maxProbability = flattenedData.Any() ? flattenedData.Max() : 0.0f;
        for (int i = 0; i < flattenedData.Count; i++)
        {
            flattenedData[i] /= maxProbability;
        }

        return flattenedData.ToArray();
    }

    public async Task UpdateDisplay()
    {
        var probabilityData = GetProbabilityData(_quantumSystem);
        Data = probabilityData;
        await JSRuntime.InvokeVoidAsync("updateThreeJsScene", probabilityData);
        StateHasChanged();
    }
}