using SkiaSharp;

namespace QuantumDiffusion3DFD.Shared;

public static class GraphicsExtensions
{
    public static float SigmoidOpacity(float probability) => (float)(1 / (1 + Math.Exp(-10 * (probability - 0.5))));

    public static string InterpolateColor(float value)
    {
        var hue = (1 - value) * 240;

        //todo play with s and l based on value as well
        var skColor = SKColor.FromHsl(hue, 100, 50);

        return $"#{skColor.Red:X2}{skColor.Green:X2}{skColor.Blue:X2}";
    }

    public static SKColor LerpColor(SKColor color1, SKColor color2, float fraction)
    {
        var r = (byte)(color1.Red + fraction * (color2.Red - color1.Red));
        var g = (byte)(color1.Green + fraction * (color2.Green - color1.Green));
        var b = (byte)(color1.Blue + fraction * (color2.Blue - color1.Blue));

        return new SKColor(r, g, b);
    }

    public static string InterpolateColorLightBlueSoftPink(float probability)
    {
        var lowColor = (R: 91, G: 206, B: 250);
        var highColor = (R: 245, G: 169, B: 184);

        var r = (byte)(lowColor.R + (highColor.R - lowColor.R) * probability);
        var g = (byte)(lowColor.G + (highColor.G - lowColor.G) * probability);
        var b = (byte)(lowColor.B + (highColor.B - lowColor.B) * probability);

        return $"#{r:X2}{g:X2}{b:X2}";
    }

    public static string InterpolateColorGreenOrange(float probability)
    {
        var greenLowColor = (R: 0, G: 255, B: 0);
        var orangeHighColor = (R: 255, G: 165, B: 0);

        var r = (byte)(greenLowColor.R + (orangeHighColor.R - greenLowColor.R) * probability);
        var g = (byte)(greenLowColor.G + (orangeHighColor.G - greenLowColor.G) * probability);
        var b = (byte)(greenLowColor.B + (orangeHighColor.B - greenLowColor.B) * probability);

        return $"#{r:X2}{g:X2}{b:X2}";
    }

    public static string InterpolateColorGrayscale(float probability)
    {
        var scaledProbability = (byte)(255 * probability);
        return $"#{scaledProbability:X2}{scaledProbability:X2}{scaledProbability:X2}";
    }

    public static string InterpolateHotCold(float probability)
    {
        var blueColdColor = (R: 0, G: 0, B: 255);
        var redHotColor = (R: 255, G: 0, B: 0);

        var r = (byte)(blueColdColor.R + (redHotColor.R - blueColdColor.R) * probability);
        var g = (byte)(blueColdColor.G + (redHotColor.G - blueColdColor.G) * probability);
        var b = (byte)(blueColdColor.B + (redHotColor.B - blueColdColor.B) * probability);

        return $"#{r:X2}{g:X2}{b:X2}";
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

        var scaledValue = probability * (colorStops.Length - 1);
        var index = (int)Math.Floor(scaledValue);
        var fraction = scaledValue - index;

        var color1 = colorStops[index];
        var color2 = colorStops[Math.Min(index + 1, colorStops.Length - 1)];

        var interpolatedColor = LerpColor(color1, color2, fraction);
        return $"#{interpolatedColor.Red:X2}{interpolatedColor.Green:X2}{interpolatedColor.Blue:X2}";
    }
}