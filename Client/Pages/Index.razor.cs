using Microsoft.AspNetCore.Components;
using Microsoft.JSInterop;
using QuantumDiffusion3DFD.Shared;

namespace QuantumDiffusion3DFD.Client.Pages;

public partial class Index
{
    private double _gaussianX0;
    private double _gaussianY0;
    private double _gaussianZ0;
    private double _gaussianSigma;
    private double _gaussianKx;
    private double _gaussianKy;
    private double _gaussianKz;

    private double _mass = 1.0;
    private double _hbar = 1.0;

    private static readonly (int x, int y, int z) Dimensions = new(10, 10, 10);
    private static float _timeStep = 0.0002f;
    private static double _spaceStep = 0.1;
    private static BoundaryType _boundaryType;
    private QuantumSystem _quantumSystem = new(Dimensions, _timeStep, _spaceStep, _boundaryType);

    private double _originalTotalEnergy;
    private double _currentTotalEnergy;
    private float[]? _currentProbabilityData;
    private float _pointSize = 21.0f;
    private int _frameCount;
    private double _frameRate;
    private DateTime _lastFrameTime = DateTime.UtcNow;

    private bool _areControlsVisible = true;
    private double _logTimeStepPosition;
    private double _logSpaceStepPosition;
    private double _hbarSliderPosition; // logarithmic scale
    private double _massSliderPosition; // logarithmic scale
    private const double MinPointLogValue = -4; // 1/10000
    private const double MaxPointLogValue = 4;  // Corresponds to 10000
    private const double MinTimeLogValue = -5; // Corresponds to 10^-5 = 0.00001
    private const double MaxTimeLogValue = -1; // Corresponds to 10^-1 = 0.1
    private const double MinSpaceLogValue = -3; // Corresponds to 10^-3 = 0.001
    private const double MaxSpaceLogValue = 1;  // Corresponds to 10^1 = 10

    private CancellationTokenSource _cancellationTokenSource = new();

    private double TimeStep => Math.Pow(10, _logTimeStepPosition);
    private double SpaceStep => Math.Pow(10, _logSpaceStepPosition);
    public static int MaxFrameRateMillis { get; set; } = 100;
    public bool Paused { get; set; } = true;
    public bool IsSimulationRunning { get; set; }

    protected override async Task OnInitializedAsync()
    {
        _gaussianX0 = Dimensions.x / 2.0;
        _gaussianY0 = Dimensions.y / 2.0;
        _gaussianZ0 = Dimensions.z / 2.0;
        _gaussianSigma = 1.0;
        _gaussianKx = 20;
        _gaussianKy = 20;
        _gaussianKz = 20;

        _cancellationTokenSource = new CancellationTokenSource();

        await ResetSimAndGraphics();
        await StartSimulation(_cancellationTokenSource.Token);
    }

    protected override async Task OnAfterRenderAsync(bool firstRender)
    {
        if (firstRender)
        {
            await JSRuntime.InvokeVoidAsync("QuantumInterop.initializeThreeJs",
                new { Dimensions.x, Dimensions.y, Dimensions.z }, _spaceStep);

            await JSRuntime.InvokeVoidAsync("QuantumInterop.updatePointSize", _pointSize);
        }
    }

    private QuantumSystem SetupQuantumSystem()
    {
        _quantumSystem = new QuantumSystem(Dimensions, _timeStep, _spaceStep, _boundaryType);

        _quantumSystem.InitializeGaussianPacket(
            _gaussianX0, _gaussianY0, _gaussianZ0, _gaussianSigma, _gaussianKx, _gaussianKy, _gaussianKz);

        _originalTotalEnergy = _quantumSystem.CalculateTotalEnergy();
        _currentTotalEnergy = _quantumSystem.CalculateTotalEnergy();

        return _quantumSystem;
    }

    private async Task ResetSimAndGraphics()
    {
        var newSystem = SetupQuantumSystem();
        _currentProbabilityData = newSystem.GetProbabilityData();
        await Update3DDisplay(_currentProbabilityData);
        StateHasChanged();
    }

    private async Task StartSimulation(CancellationToken token)
    {
        IsSimulationRunning = true;

        while (!token.IsCancellationRequested)
        {
            if (Paused) await Task.Delay(100, token); // Small delay to reduce CPU usage
            else
            {
                await Task.Run(async () =>
                {
                    _quantumSystem.ApplySingleTimeEvolutionStep();
                    _currentTotalEnergy = _quantumSystem.CalculateTotalEnergy();
                    _currentProbabilityData = _quantumSystem.GetProbabilityData();
                    await Update3DDisplay(_currentProbabilityData);
                }, token);

                if (token.IsCancellationRequested) break;
                await DelayUntilNextFrame(token);
                await InvokeAsync(StateHasChanged);
            }
        }

        IsSimulationRunning = false;
    }

    private async Task DelayUntilNextFrame(CancellationToken token)
    {
        var currentFrameTime = DateTime.UtcNow;
        var elapsed = currentFrameTime - _lastFrameTime;

        if (elapsed.TotalSeconds > 0)
            _frameRate = 1.0 / elapsed.TotalSeconds;
        _lastFrameTime = currentFrameTime;
        _frameCount++;

        var difference = MaxFrameRateMillis - (int)elapsed.TotalMilliseconds;
        if (difference > 0)
            await Task.Delay(difference, token);
    }

    public async Task Update3DDisplay(float[]? probabilityData)
    {
        probabilityData ??= _quantumSystem.GetProbabilityData();
        await JSRuntime.InvokeVoidAsync("QuantumInterop.updateThreeJsScene", probabilityData);
    }

    private async Task RestartSimulation()
    {
        if (IsSimulationRunning)
        {
            _cancellationTokenSource.Cancel();
            _cancellationTokenSource.Dispose();
            _cancellationTokenSource = new CancellationTokenSource();
        }

        await ResetSimAndGraphics();
        await StartSimulation(_cancellationTokenSource.Token);
    }

    private void TogglePause() => Paused = !Paused;

    private void ToggleControlsVisibility() => _areControlsVisible = !_areControlsVisible;

    private string GetControlPanelClass() => _areControlsVisible ? "expanded" : "collapsed";

    private async Task UpdateProbabilitySphereScale(ChangeEventArgs e) => 
        await JSRuntime.InvokeVoidAsync("QuantumInterop.updatePointSize", 
            float.Parse(e.Value?.ToString() ?? "1"));

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
}