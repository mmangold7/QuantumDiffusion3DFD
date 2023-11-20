using Microsoft.AspNetCore.Components;
using Microsoft.JSInterop;
using QuantumDiffusion3DFD.Shared;
using SkiaSharp;
using System.Diagnostics;

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
    private int _frameCount;
    private float _frameRate;
    private DateTime _lastFrameTime = DateTime.UtcNow;
    private float[]? _previousProbabilityData;
    private float[]? _currentProbabilityData;
    private QuantumSystem? _quantumSystem;
    private CancellationTokenSource _cancellationTokenSource = new();
    private Stopwatch? _stopwatch;

    private bool Paused { get; set; }
    private int MaxFramesPerSecond { get; set; } = 30;
    private float TimeStep { get; set; } = 0.0001f;
    private float GaussianX0 { get; set; } = BoxX / 2.0f;
    private float GaussianY0 { get; set; } = BoxY / 2.0f;
    private float GaussianZ0 { get; set; } = BoxZ / 2.0f;
    private float GaussianSigma { get; set; } = 1;
    private float GaussianKx { get; set; } = 20;
    private float GaussianKy { get; set; } = -50;
    private float GaussianKz { get; set; } = 75;
    private float Mass { get; set; } = 10;
    private float Hbar { get; set; } = 10;
    private BoundaryType BoundaryType { get; set; } = BoundaryType.Reflective;
    private float OriginalTotalEnergy { get; set; }
    private float CurrentTotalEnergy { get; set; }

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
        Console.WriteLine($"User {(success ? "accepted" : "dismissed")} the A2HS prompt");
    }

    private QuantumSystem FromUiValues() =>
        new(new ValueTuple<int, int, int>(BoxX, BoxY, BoxZ), BoundaryType, TimeStep, SpaceStep, Mass, Hbar);

    private async Task ResetSimAndGraphics()
    {
        _frameCount = 0;
        _frameRate = 0;

        OriginalTotalEnergy = 0;
        CurrentTotalEnergy = 0;

        _quantumSystem = FromUiValues();
        _quantumSystem.InitializeGaussianPacket(GaussianX0, GaussianY0, GaussianZ0, GaussianSigma, GaussianKx, GaussianKy, GaussianKz);

        _previousProbabilityData = null;
        _currentProbabilityData = _quantumSystem.GetNormalizedProbabilityData();

        await Update3DDisplay(_currentProbabilityData);

        StateHasChanged();
    }

    private async Task StartSimulation(CancellationToken token)
    {
        _isSimulationRunning = true;

        if (Extensions.IsDebug) _stopwatch = Stopwatch.StartNew();

        while (!token.IsCancellationRequested)
        {
            if (Paused) await Task.Delay(100, token); // Small delay to reduce CPU usage
            else
            {
                await Task.Run(async () =>
                {
                    _quantumSystem.ApplySingleTimeEvolutionStepEuler();
                    if (Extensions.IsDebug) Extensions.LogMethodTime(nameof(_quantumSystem.ApplySingleTimeEvolutionStepEuler), _stopwatch);
                    
                    CurrentTotalEnergy = _quantumSystem.CalculateTotalEnergy();
                    if (OriginalTotalEnergy == 0) OriginalTotalEnergy = CurrentTotalEnergy;
                    //if (Extensions.IsDebug) Extensions.LogMethodTime(nameof(_quantumSystem.CalculateTotalEnergy), _stopwatch);
                    
                    _currentProbabilityData = _quantumSystem.GetNormalizedProbabilityData();
                    //if (Extensions.IsDebug) Extensions.LogMethodTime(nameof(_quantumSystem.GetNormalizedProbabilityData), _stopwatch);
                    
                    await Update3DDisplay(_currentProbabilityData);
                    if (Extensions.IsDebug) Extensions.LogMethodTime(nameof(Update3DDisplay), _stopwatch);
                }, token);

                if (token.IsCancellationRequested) break;
                
                await DelayUntilNextFrame(token);
                //if (Extensions.IsDebug) Extensions.LogMethodTime(nameof(DelayUntilNextFrame), _stopwatch);
                
                await InvokeAsync(StateHasChanged);
                //if (Extensions.IsDebug) Extensions.LogMethodTime(nameof(StateHasChanged), _stopwatch);
            }
        }

        _isSimulationRunning = false;
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

    public async Task Update3DDisplay(float[] probabilityData)
    {
        var updatedData = new List<object>();
        float maxProbability = 1;
        float updateThreshold = maxProbability * 0.1f;
        float opacityScale = 0.75f;

        for (int i = 0; i < probabilityData.Length; i++)
        {
            var newProbability = probabilityData[i];
            if (_previousProbabilityData == null || Math.Abs(newProbability - _previousProbabilityData[i]) > updateThreshold)
            {
                var color = InterpolateColor(newProbability);
                var opacity = SigmoidOpacity(newProbability * opacityScale);
                updatedData.Add(new { index = i, color, opacity });
                if (_previousProbabilityData != null) _previousProbabilityData[i] = newProbability;
            }
        }

        if (updatedData.Count > 0)
        {
            await JSRuntime.InvokeVoidAsync("QuantumInterop.updateThreeJsScene", updatedData);
        }
        _previousProbabilityData ??= probabilityData;
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
        if (_isSimulationRunning)
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
}