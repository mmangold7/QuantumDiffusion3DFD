using Microsoft.AspNetCore.Components;
using Microsoft.JSInterop;
using QuantumDiffusion3DFD.Shared;
using SkiaSharp;
using System.Diagnostics;

namespace QuantumDiffusion3DFD.Client.Pages;

public partial class Index
{
    // shell command for exposing api to network: "npx iisexpress-proxy https://localhost:7223 to 3000"

    private double _gaussianX0;
    private double _gaussianY0;
    private double _gaussianZ0;
    private double _gaussianSigma;
    private double _gaussianKx;
    private double _gaussianKy;
    private double _gaussianKz;

    private const int cubeDimensionSize = 10;

    private static int _boxX = cubeDimensionSize;
    private static int _boxY = cubeDimensionSize;
    private static int _boxZ = cubeDimensionSize;

    private static (int x, int y, int z) Dimensions => new(_boxX, _boxY, _boxZ);
    private static double _mass = 1.0;
    private static double _hbar = 2.0;
    private static float _timeStep = 0.0001f;
    private static double _spaceStep = 0.1;
    private static BoundaryType _boundaryType;
    private QuantumSystem _quantumSystem = BuildSystemWithUiParameters();

    private static QuantumSystem BuildSystemWithUiParameters() =>
        new(Dimensions, _boundaryType, _timeStep, _spaceStep, _mass, _hbar);

    private double _originalTotalEnergy;
    private double _currentTotalEnergy;
    private float[]? previousProbabilityData;
    private float[]? _currentProbabilityData;
    private float _pointSize = 21.0f;
    private int _frameCount;
    private double _frameRate;
    private DateTime _lastFrameTime = DateTime.UtcNow;
    private static int _maxFramesPerSecond = 20;

    private bool _areControlsVisible = true;
    private double _logTimeStepPosition;
    private double _logSpaceStepPosition;
    private double _hbarSliderPosition; // logarithmic scale
    private double _massSliderPosition; // logarithmic scale
    private const double MinPointLogValue = -4; // 1/10000
    private const double MaxPointLogValue = 4;  // Corresponds to 10000
    //private const double MinTimeLogValue = -5; // Corresponds to 10^-5 = 0.00001
    //private const double MaxTimeLogValue = -1; // Corresponds to 10^-1 = 0.1
    //private const double MinSpaceLogValue = -3; // Corresponds to 10^-3 = 0.001
    //private const double MaxSpaceLogValue = 1;  // Corresponds to 10^1 = 10

    private CancellationTokenSource _cancellationTokenSource = new();
    Stopwatch? stopwatch;

    private double TimeStep => Math.Pow(10, _logTimeStepPosition);
    private double SpaceStep => Math.Pow(10, _logSpaceStepPosition);
    public static int MaxFrameRateMillis => 1000 / _maxFramesPerSecond;
    public bool Paused { get; set; } = true;
    public bool IsSimulationRunning { get; set; }

    protected override async Task OnInitializedAsync()
    {
        _hbarSliderPosition = _hbar;
        _massSliderPosition = _mass;
        _gaussianX0 = Dimensions.x / 2.0;
        _gaussianY0 = Dimensions.y / 2.0;
        _gaussianZ0 = Dimensions.z / 2.0;
        _gaussianSigma = .75;
        _gaussianKx = 10;
        _gaussianKy = 10;
        _gaussianKz = 10;
    }

    protected override async Task OnAfterRenderAsync(bool firstRender)
    {
        if (firstRender)
        {
            await JSRuntime.InvokeVoidAsync("QuantumInterop.initializeThreeJs",
                new { Dimensions.x, Dimensions.y, Dimensions.z }, _spaceStep);

            await JSRuntime.InvokeVoidAsync("QuantumInterop.updatePointSize", _pointSize);

            await RestartSimulation();
        }
    }

    private async Task InstallApp()
    {
        var success = await JSRuntime.InvokeAsync<bool>("showPWAInstallPrompt");
        Console.WriteLine($"User {(success ? "accepted" : "dismissed")} the A2HS prompt");
    }

    private QuantumSystem SetupQuantumSystem()
    {
        _quantumSystem = BuildSystemWithUiParameters();
        _quantumSystem.InitializeGaussianPacket(_gaussianX0, _gaussianY0, _gaussianZ0, _gaussianSigma, _gaussianKx, _gaussianKy, _gaussianKz);
        return _quantumSystem;
    }

    private async Task ResetSimAndGraphics()
    {
        var newSystem = SetupQuantumSystem();
        _currentProbabilityData = newSystem.GetNormalizedProbabilityData();
        _frameCount = 0;
        _frameRate = 0;
        _originalTotalEnergy = 0;
        _currentTotalEnergy = 0;
        await Update3DDisplay(_currentProbabilityData);
        StateHasChanged();
    }

    private async Task StartSimulation(CancellationToken token)
    {
        IsSimulationRunning = true;

        var debugOutput = false;
        if (debugOutput) stopwatch = Stopwatch.StartNew();

        while (!token.IsCancellationRequested)
        {
            if (Paused) await Task.Delay(100, token); // Small delay to reduce CPU usage
            else
            {
                await Task.Run(async () =>
                {
                    _quantumSystem.ApplySingleTimeEvolutionStepEuler();
                    if (debugOutput) Extensions.LogMethodTime(nameof(_quantumSystem.ApplySingleTimeEvolutionStepEuler), stopwatch);
                    
                    _currentTotalEnergy = _quantumSystem.CalculateTotalEnergy();
                    if (_originalTotalEnergy == 0) _originalTotalEnergy = _currentTotalEnergy;
                    if (debugOutput) Extensions.LogMethodTime(nameof(_quantumSystem.CalculateTotalEnergy), stopwatch);
                    
                    _currentProbabilityData = _quantumSystem.GetNormalizedProbabilityData();
                    if (debugOutput) Extensions.LogMethodTime(nameof(_quantumSystem.GetNormalizedProbabilityData), stopwatch);
                    
                    await Update3DDisplay(_currentProbabilityData);
                    if (debugOutput) Extensions.LogMethodTime(nameof(Update3DDisplay), stopwatch);
                }, token);

                if (token.IsCancellationRequested) break;
                
                await DelayUntilNextFrame(token);
                if (debugOutput) Extensions.LogMethodTime(nameof(DelayUntilNextFrame), stopwatch);
                
                await InvokeAsync(StateHasChanged);
                if (debugOutput) Extensions.LogMethodTime(nameof(StateHasChanged), stopwatch);
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

    public async Task Update3DDisplay(float[] probabilityData)
    {
        var updatedData = new List<object>();
        float maxProbability = 1;
        float updateThreshold = maxProbability * 0.1f;
        float opacityScale = 0.75f;

        for (int i = 0; i < probabilityData.Length; i++)
        {
            var newProbability = probabilityData[i];
            if (previousProbabilityData == null || Math.Abs(newProbability - previousProbabilityData[i]) > updateThreshold)
            {
                var color = InterpolateColorGreenOrange(newProbability);
                var opacity = SigmoidOpacity(newProbability * opacityScale);
                updatedData.Add(new { index = i, color, opacity });
                if (previousProbabilityData != null) previousProbabilityData[i] = newProbability;
            }
        }

        if (updatedData.Count > 0)
        {
            await JSRuntime.InvokeVoidAsync("QuantumInterop.updateThreeJsScene", updatedData);
        }
        previousProbabilityData ??= probabilityData;
    }

    private string InterpolateColorLightBlueSoftPink(float probability)
    {
        // Starting color: Light Blue (#5BCEFA)
        var lowColor = (R: 91, G: 206, B: 250);
        // Ending color: Soft Pink (#F5A9B8)
        var highColor = (R: 245, G: 169, B: 184);

        // Linear interpolation between the colors based on probability
        byte r = (byte)(lowColor.R + (highColor.R - lowColor.R) * probability);
        byte g = (byte)(lowColor.G + (highColor.G - lowColor.G) * probability);
        byte b = (byte)(lowColor.B + (highColor.B - lowColor.B) * probability);

        return $"#{r:X2}{g:X2}{b:X2}";
    }

    private string InterpolateColorGreenOrange(float probability)
    {
        // Define the endpoint colors: Green (low) and Orange (high)
        var lowColor = (R: 0, G: 255, B: 0); // Green
        var highColor = (R: 255, G: 165, B: 0); // Orange

        // Linear interpolation between the colors based on probability
        byte r = (byte)(lowColor.R + (highColor.R - lowColor.R) * probability);
        byte g = (byte)(lowColor.G + (highColor.G - lowColor.G) * probability);
        byte b = (byte)(lowColor.B + (highColor.B - lowColor.B) * probability);

        return $"#{r:X2}{g:X2}{b:X2}";
    }


    private string InterpolateColorGrayscale(float probability)
    {
        var scaledProbability = (byte)(255 * probability);
        return $"#{scaledProbability:X2}{scaledProbability:X2}{scaledProbability:X2}";
    }

    private string InterpolateHotCold(float probability)
    {
        // Endpoint colors
        var coldColor = (R: 0, G: 0, B: 255); // Blue
        var hotColor = (R: 255, G: 0, B: 0); // Red

        // Interpolate
        byte r = (byte)(coldColor.R + (hotColor.R - coldColor.R) * probability);
        byte g = (byte)(coldColor.G + (hotColor.G - coldColor.G) * probability);
        byte b = (byte)(coldColor.B + (hotColor.B - coldColor.B) * probability);

        return $"#{r:X2}{g:X2}{b:X2}";
    }

    private string InterpolateColor(float value)
    {
        // Convert value to HSL color (hue from blue to red)
        float hue = (1 - value) * 240;

        // Using SkiaSharp to convert HSL to RGB
        var skColor = SKColor.FromHsl(hue, 100, 50); // S and L values are percentages

        // Format as hex string
        return $"#{skColor.Red:X2}{skColor.Green:X2}{skColor.Blue:X2}";
    }

    private SKColor LerpColor(SKColor color1, SKColor color2, float fraction)
    {
        byte r = (byte)(color1.Red + fraction * (color2.Red - color1.Red));
        byte g = (byte)(color1.Green + fraction * (color2.Green - color1.Green));
        byte b = (byte)(color1.Blue + fraction * (color2.Blue - color1.Blue));

        return new SKColor(r, g, b);
    }

    private string InterpolateColorRainbow(float probability)
    {
        var colorStops = new[]
        {
            SKColors.Red,
            SKColors.Orange,
            SKColors.Yellow,
            SKColors.Green,
            SKColors.Blue,
            SKColors.Indigo,
            SKColors.Violet
        };

        float scaledValue = probability * (colorStops.Length - 1);
        int index = (int)Math.Floor(scaledValue);
        float frac = scaledValue - index;

        var color1 = colorStops[index];
        var color2 = colorStops[Math.Min(index + 1, colorStops.Length - 1)];

        var interpolatedColor = LerpColor(color1, color2, frac);
        return $"#{interpolatedColor.Red:X2}{interpolatedColor.Green:X2}{interpolatedColor.Blue:X2}";
    }

    private float SigmoidOpacity(float probability) => (float)(1 / (1 + Math.Exp(-10 * (probability - 0.5))));

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
    }

    private void OnMassSliderChanged(ChangeEventArgs e)
    {
        _massSliderPosition = Convert.ToDouble(e.Value);
        _mass = Math.Pow(10, _massSliderPosition);
    }
}