using System.Diagnostics;

namespace QuantumDiffusion3DFD.Shared;

public static class Extensions
{
    // shell command for exposing api to network: "npx iisexpress-proxy https://localhost:7223 to 3000"

    public static void LogMethodTime(string methodName, Stopwatch? stopwatch)
    {
        Console.WriteLine($"{methodName} took: {stopwatch?.ElapsedMilliseconds} ms");
        stopwatch?.Restart();
    }

#if DEBUG
    public const bool IsDebug = false;
#else
    public const bool IsDebug = false;
#endif
}