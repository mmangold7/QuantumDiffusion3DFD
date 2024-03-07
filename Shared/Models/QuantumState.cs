namespace QuantumDiffusion3DFD.Shared.Models;

public class QuantumState
{
    public float OriginalTotalEnergy { get; set; }
    public float CurrentTotalEnergy { get; set; }
    public List<object> ProbabilityData { get; set; }
}