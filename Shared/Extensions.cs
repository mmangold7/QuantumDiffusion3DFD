using System.Diagnostics;

namespace QuantumDiffusion3DFD.Shared;

public static class Extensions
{
    public static void LogMethodTime(string methodName, Stopwatch? stopwatch)
    {
        Console.WriteLine($"{methodName} took: {stopwatch?.ElapsedMilliseconds} ms");
        stopwatch?.Restart();
    }
}