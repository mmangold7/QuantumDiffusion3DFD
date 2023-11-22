namespace QuantumDiffusion3DFD.Shared;

public class QuantumState
{
    public float OriginalTotalEnergy { get; set; }
    public float CurrentTotalEnergy { get; set; }
    public float[] ProbabilityData { get; set; }
    public List<object> SignificantlyChangedProbabilityData { get; set; }
}