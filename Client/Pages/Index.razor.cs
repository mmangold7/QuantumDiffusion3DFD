using QuantumDiffusion3DFD.Shared;
using System.Numerics;
using Microsoft.JSInterop;
using Microsoft.AspNetCore.Components;

namespace QuantumDiffusion3DFD.Client.Pages;

public partial class Index
{
    private QuantumSystem _quantumSystem;
    private Complex _wavefunctionValue;
    private int _frameCount;
    private double _frameRate;
    private (int x, int y, int z) dimensions = new(10, 10, 10);
    private double spacing = 0.1;
    private double _hbar = 1.0;
    private double _mass = 1.0;
    private double _hbarSliderPosition = 0; // logarithmic scale
    private double _massSliderPosition = 0; // logarithmic scale
    private const double MinLogValue = -4; // Corresponds to 1/10000
    private const double MaxLogValue = 4;  // Corresponds to 10000
    private float _pointSize = 13.0f;

    public bool Paused { get; set; }
    public bool RunningSimulation { get; set; }
    public static int SimulationDelayMilliseconds { get; set; } = 33;

    private async Task UpdatePointSize(ChangeEventArgs e)
    {
        _pointSize = float.Parse(e.Value.ToString());
        await JSRuntime.InvokeVoidAsync("updatePointSize", _pointSize);
    }

    private void OnHbarSliderChanged(ChangeEventArgs e)
    {
        _hbarSliderPosition = Convert.ToDouble(e.Value);
        _hbar = Math.Pow(10, _hbarSliderPosition);
        _quantumSystem.Hbar = _hbar;
    }

    private void OnMassSliderChanged(ChangeEventArgs e)
    {
        _massSliderPosition = Convert.ToDouble(e.Value);
        _mass = Math.Pow(10, _massSliderPosition);
        _quantumSystem.SingleParticleMass = _mass;
    }

    private void TogglePause() => Paused = !Paused;

    protected override async Task OnAfterRenderAsync(bool firstRender)
    {
        if (firstRender)
        {
            await JSRuntime.InvokeVoidAsync("initializeThreeJs", new { x = dimensions.x, y = dimensions.y, z = dimensions.z }, spacing);
            await JSRuntime.InvokeVoidAsync("updatePointSize", _pointSize);
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
        _quantumSystem = new QuantumSystem(dimensions, 0.0002, spacing, BoundaryType.Reflective);
        _quantumSystem.InitializeGaussianPacket(_quantumSystem.Dimensions.x / 2, _quantumSystem.Dimensions.y / 2, _quantumSystem.Dimensions.z / 2, 1, 0, 0, 0);
    }

    private async Task StartQuantumSimulationLoop(CancellationToken token, bool paused = true)
    {
        RunningSimulation = true;
        Paused = paused;

        DateTime lastFrameTime = DateTime.UtcNow;

        while (RunningSimulation)
        {
            DateTime currentFrameTime = DateTime.UtcNow;
            var elapsed = (currentFrameTime - lastFrameTime);

            if (!Paused)
            {
                _frameRate = 1.0 / elapsed.TotalSeconds;
                lastFrameTime = currentFrameTime;
                _frameCount++;

                _quantumSystem.ApplySingleTimeEvolutionStep();

                UpdateTextDisplay();
                await Update3DDisplay();

                StateHasChanged();
            }

            var difference = SimulationDelayMilliseconds - elapsed.Milliseconds;
            if (difference > 0)
                await Task.Delay(difference, token);
            else
                await Task.Delay(SimulationDelayMilliseconds, token);
        }
    }

    private void UpdateTextDisplay()
    {
        _wavefunctionValue = _quantumSystem.GetWavefunctionValue(dimensions.x / 2, dimensions.y / 2, dimensions.z / 2);
    }

    public float[] GetProbabilityData(QuantumSystem quantumSystem)
    {
        var probabilityDensity = quantumSystem.CalculateProbabilityDensity();
        var flattenedData = new List<float>();

        for (int x = 0; x < _quantumSystem.Dimensions.x; x++)
        {
            for (int y = 0; y < _quantumSystem.Dimensions.y; y++)
            {
                for (int z = 0; z < _quantumSystem.Dimensions.z; z++)
                {
                    float probability = (float)probabilityDensity[x, y, z];

                    // Check for invalid numbers
                    if (float.IsNaN(probability) || float.IsInfinity(probability))
                    {
                        Logger.LogWarning($"Invalid probability value detected at ({x}, {y}, {z}): {probability}");
                        probability = 0;  // Or handle this case as appropriate
                    }

                    flattenedData.Add(probability);
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

    public async Task Update3DDisplay()
    {
        var probabilityData = GetProbabilityData(_quantumSystem);
        await JSRuntime.InvokeVoidAsync("updateThreeJsScene", probabilityData);
    }
}