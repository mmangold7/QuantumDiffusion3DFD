using SkiaSharp;

namespace QuantumDiffusion3DFD.Shared;

public static class GraphicsExtensions
{
    public static string InterpolateColorLightBlueSoftPink(float probability)
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

    public static string InterpolateColorGreenOrange(float probability)
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


    public static string InterpolateColorGrayscale(float probability)
    {
        var scaledProbability = (byte)(255 * probability);
        return $"#{scaledProbability:X2}{scaledProbability:X2}{scaledProbability:X2}";
    }

    public static string InterpolateHotCold(float probability)
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

    public static string InterpolateColor(float value)
    {
        // Convert value to HSL color (hue from blue to red)
        float hue = (1 - value) * 240;

        // Using SkiaSharp to convert HSL to RGB
        var skColor = SKColor.FromHsl(hue, 100, 50); // S and L values are percentages

        // Format as hex string
        return $"#{skColor.Red:X2}{skColor.Green:X2}{skColor.Blue:X2}";
    }

    public static SKColor LerpColor(SKColor color1, SKColor color2, float fraction)
    {
        byte r = (byte)(color1.Red + fraction * (color2.Red - color1.Red));
        byte g = (byte)(color1.Green + fraction * (color2.Green - color1.Green));
        byte b = (byte)(color1.Blue + fraction * (color2.Blue - color1.Blue));

        return new SKColor(r, g, b);
    }

    public static string InterpolateColorRainbow(float probability)
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

    public static float SigmoidOpacity(float probability) => (float)(1 / (1 + Math.Exp(-10 * (probability - 0.5))));
}