using Microsoft.AspNetCore.Components;
using Microsoft.JSInterop;
using QuantumDiffusion3DFD.Shared;

namespace QuantumDiffusion3DFD.Client.Pages;

public partial class Index
{
    private QuantumSystem _quantumSystem;
    private int _frameCount;
    private double _frameRate;
    private (int x, int y, int z) dimensions = new(10, 10, 10);
    private double _hbar = 1.0;
    private double _mass = 1.0;
    private double _hbarSliderPosition = 0; // logarithmic scale
    private double _massSliderPosition = 0; // logarithmic scale
    private const double MinLogValue = -4; // Corresponds to 1/10000
    private const double MaxLogValue = 4;  // Corresponds to 10000
    private float _pointSize = 13.0f;
    private CancellationToken simulationToken;
    private CancellationTokenSource? simulationTokenSource;
    private float _timeStep = 0.0002f;
    private double _spaceStep = 0.1;

    private double _logTimeStepPosition;
    private double _logSpaceStepPosition;
    private const double MinTimeLogValue = -5; // Corresponds to 10^-5 = 0.00001
    private const double MaxTimeLogValue = -1; // Corresponds to 10^-1 = 0.1
    private const double MinSpaceLogValue = -3; // Corresponds to 10^-3 = 0.001
    private const double MaxSpaceLogValue = 1;  // Corresponds to 10^1 = 10

    private double _originalTotalEnergy;
    private double _currentTotalEnergy;

    private double TimeStep => Math.Pow(10, _logTimeStepPosition);
    private double SpaceStep => Math.Pow(10, _logSpaceStepPosition);

    private BoundaryType _boundaryType;

    private double _gaussianX0;
    private double _gaussianY0;
    private double _gaussianZ0 ;
    private double _gaussianSigma = 1.0; // Default value
    private double _gaussianKx;
    private double _gaussianKy;
    private double _gaussianKz;

    public bool Paused { get; set; } = true;
    public bool RunningSimulation { get; set; }
    public static int SimulationDelayMilliseconds { get; set; } = 33;

    private async Task UpdatePointSize(ChangeEventArgs e)
    {
        _pointSize = float.Parse(e.Value.ToString());
        await JSRuntime.InvokeVoidAsync("updatePointSize", _pointSize);
    }

    private void OnTimeStepChanged(ChangeEventArgs e)
    {
        _logTimeStepPosition = Convert.ToDouble(e.Value);
        _timeStep = (float)TimeStep;
    }

    private void OnSpaceStepChanged(ChangeEventArgs e)
    {
        _logSpaceStepPosition = Convert.ToDouble(e.Value);
        _spaceStep = (float)SpaceStep;
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

    private async Task RestartSimulation()
    {
        if (simulationToken != null)
            simulationTokenSource?.Cancel();

        simulationTokenSource = new CancellationTokenSource();
        simulationToken = simulationTokenSource.Token;
        
        SetupQuantumSystem();
        _quantumSystem.InitializeGaussianPacket(_gaussianX0, _gaussianY0, _gaussianZ0, _gaussianSigma, _gaussianKx, _gaussianKy, _gaussianKz);
        StateHasChanged();

        await StartQuantumSimulationLoop(simulationToken, Paused);
    }

    private void TogglePause() => Paused = !Paused;

    protected override async Task OnAfterRenderAsync(bool firstRender)
    {
        if (firstRender)
        {
            await JSRuntime.InvokeVoidAsync("initializeThreeJs", new { x = dimensions.x, y = dimensions.y, z = dimensions.z }, _spaceStep);
            await JSRuntime.InvokeVoidAsync("updatePointSize", _pointSize);
        }
    }

    protected override async Task OnInitializedAsync()
    {
        await RestartSimulation();
    }

    private void SetupQuantumSystem()
    {
        _quantumSystem = new QuantumSystem(dimensions, _timeStep, _spaceStep, _boundaryType);

        _gaussianX0 = _quantumSystem.Dimensions.x / 2;
        _gaussianY0 = _quantumSystem.Dimensions.y / 2;
        _gaussianZ0 = _quantumSystem.Dimensions.z / 2;
        _gaussianSigma = 1.0; // Default value
        _gaussianKx = 10;
        _gaussianKy = 10;
        _gaussianKz = 10;

        _quantumSystem.InitializeGaussianPacket(_gaussianX0, _gaussianY0, _gaussianZ0, _gaussianSigma, _gaussianKx, _gaussianKy, _gaussianKz);
        _originalTotalEnergy = _quantumSystem.CalculateTotalEnergy();
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
                _currentTotalEnergy = _quantumSystem.CalculateTotalEnergy();

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

    //private async Task ToggleControlPanel()
    //{
    //    await JSRuntime.InvokeVoidAsync("toggleControls");
    //}

    private bool _areControlsVisible = true; // Controls are hidden by default

    private string ToggleText => _areControlsVisible ? "Hide Controls" : "Show Controls";

    private void ToggleControlsVisibility()
    {
        _areControlsVisible = !_areControlsVisible;
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