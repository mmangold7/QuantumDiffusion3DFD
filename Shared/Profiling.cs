using System.Diagnostics;

namespace QuantumDiffusion3DFD.Shared;

public static class Profiling
{
    // shell command for exposing api to network: "npx iisexpress-proxy https://localhost:7223 to 3000"

    public const bool ShouldLogMethodProfiles = true;

    public const string DebugEnvironmentName = "DEVELOPMENT";
    public const string ReleaseEnvironmentName = "PRODUCTION";

#if RELEASE
    public static string EnvironmentName => ReleaseEnvironmentName;
#else
    public static string EnvironmentName => DebugEnvironmentName;
#endif

    public static bool IsDebug = EnvironmentName == DebugEnvironmentName;
    public static bool IsRelease = EnvironmentName == ReleaseEnvironmentName;

    public static void RunWithClockingLog(this Action action)
    {
        if (ShouldLogMethodProfiles && IsDebug)
        {
            var methodName = GetCallingMethodName();
            var stopwatch = Stopwatch.StartNew();
            action();
            stopwatch.Stop();
            LogMethodTime(methodName, stopwatch.ElapsedMilliseconds);
        }
        else
        {
            action();
        }
    }

    public static async Task RunWithClockingLogAsync(this Func<Task> action)
    {
        if (ShouldLogMethodProfiles && IsDebug)
        {
            var methodName = GetCallingMethodName();
            var stopwatch = Stopwatch.StartNew();
            await action();
            stopwatch.Stop();
            LogMethodTime(methodName, stopwatch.ElapsedMilliseconds);
        }
        else
        {
            await action();
        }
    }

    private static string GetCallingMethodName()
    {
        var frame = new StackFrame(2);
        var method = frame.GetMethod();
        return method.Name;
    }

    private static void LogMethodTime(string methodName, long elapsedMilliseconds) => 
        Console.WriteLine($"{methodName} took: {elapsedMilliseconds} ms");
}