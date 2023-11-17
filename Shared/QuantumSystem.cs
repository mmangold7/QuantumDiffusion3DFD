using System.Numerics;

namespace QuantumDiffusion3DFD.Shared;

public class QuantumSystem
{
    private Complex[,,] _wavefunction;
    private readonly Complex[,,] _laplacianWavefunction;
    private readonly double[,,] _potential;
    private readonly double _timeStep;
    private readonly double _deltaX;
    private readonly BoundaryType _boundaryType;

    public double Hbar { get; set; } = 1.0;
    public double SingleParticleMass { get; set; } = 1.0;
    public (int x, int y, int z) Dimensions { get; set; }

    public QuantumSystem((int x, int y, int z) dimensions, double timeStep, double spaceStep, BoundaryType boundaryType)
    {
        Dimensions = dimensions;

        _timeStep = timeStep;
        _deltaX = spaceStep;
        _boundaryType = boundaryType;

        _wavefunction = new Complex[dimensions.x, dimensions.y, dimensions.z];
        _laplacianWavefunction = new Complex[dimensions.x, dimensions.y, dimensions.z];
        _potential = new double[dimensions.x, dimensions.y, dimensions.z];
    }

    public double CalculateTotalEnergy()
    {
        var totalEnergy = 0.0;

        for (var x = 0; x < Dimensions.x; x++)
        {
            for (var y = 0; y < Dimensions.y; y++)
            {
                for (var z = 0; z < Dimensions.z; z++)
                {
                    var psi = _wavefunction[x, y, z];
                    var probabilityDensity = psi.Magnitude * psi.Magnitude;

                    var laplacianPsi = _laplacianWavefunction[x, y, z]; // Use the stored Laplacian
                    var kineticEnergy = -(Hbar * Hbar / (2 * SingleParticleMass)) * (laplacianPsi * Complex.Conjugate(psi)).Real;

                    var potentialEnergy = _potential[x, y, z] * probabilityDensity;
                    totalEnergy += kineticEnergy + potentialEnergy;
                }
            }
        }

        return totalEnergy;
    }

    public void ApplySingleTimeEvolutionStep()
    {
        var newWavefunction = new Complex[Dimensions.x, Dimensions.y, Dimensions.z];

        for (var x = 0; x < Dimensions.x; x++)
        {
            for (var y = 0; y < Dimensions.y; y++)
            {
                for (var z = 0; z < Dimensions.z; z++)
                {
                    _laplacianWavefunction[x, y, z] = CalculateLaplacian(_wavefunction, x, y, z, _boundaryType);

                    var timeDerivative = (-Hbar * Hbar / (2 * SingleParticleMass)) * _laplacianWavefunction[x, y, z] + _potential[x, y, z] * _wavefunction[x, y, z];
                    newWavefunction[x, y, z] = _wavefunction[x, y, z] - (Complex.ImaginaryOne / Hbar) * timeDerivative * _timeStep;
                }
            }
        }

        _wavefunction = newWavefunction;
    }

    private Complex CalculateLaplacian(Complex[,,] wavefunction, int x, int y, int z, BoundaryType boundaryType)
    {
        var dx2 = _deltaX * _deltaX; // Assuming deltaX is the grid spacing

        var d2PsiDx2 = GetSecondDerivative(wavefunction, x, y, z, 0, boundaryType, dx2);
        var d2PsiDy2 = GetSecondDerivative(wavefunction, x, y, z, 1, boundaryType, dx2);
        var d2PsiDz2 = GetSecondDerivative(wavefunction, x, y, z, 2, boundaryType, dx2);

        return d2PsiDx2 + d2PsiDy2 + d2PsiDz2;
    }

    private Complex GetSecondDerivative(Complex[,,] wavefunction, int x, int y, int z, int dimension, BoundaryType boundaryType, double dx2)
    {
        int plusIndex, minusIndex;
        Complex valuePlus, valueMinus;

        switch (dimension)
        {
            case 0: // x-dimension
                plusIndex = Mod(x + 1, Dimensions.x, boundaryType);
                minusIndex = Mod(x - 1, Dimensions.x, boundaryType);
                valuePlus = plusIndex != -1 ? wavefunction[plusIndex, y, z] : 0;
                valueMinus = minusIndex != -1 ? wavefunction[minusIndex, y, z] : 0;
                break;
            case 1: // y-dimension
                plusIndex = Mod(y + 1, Dimensions.y, boundaryType);
                minusIndex = Mod(y - 1, Dimensions.y, boundaryType);
                valuePlus = plusIndex != -1 ? wavefunction[x, plusIndex, z] : 0;
                valueMinus = minusIndex != -1 ? wavefunction[x, minusIndex, z] : 0;
                break;
            case 2: // z-dimension
                plusIndex = Mod(z + 1, Dimensions.z, boundaryType);
                minusIndex = Mod(z - 1, Dimensions.z, boundaryType);
                valuePlus = plusIndex != -1 ? wavefunction[x, y, plusIndex] : 0;
                valueMinus = minusIndex != -1 ? wavefunction[x, y, minusIndex] : 0;
                break;
            default:
                throw new ArgumentException("Invalid dimension");
        }

        return (valuePlus - 2 * wavefunction[x, y, z] + valueMinus) / dx2;
    }

    private static int Mod(int a, int b, BoundaryType boundaryType)
    {
        switch (boundaryType)
        {
            case BoundaryType.Reflective:
                if (a < 0 || a >= b)
                    return Math.Abs(b - Math.Abs(a) % b) % b;
                break;
            case BoundaryType.Absorbing:
                if (a < 0 || a >= b)
                    return -1;
                break;
            default:
                throw new ArgumentOutOfRangeException(nameof(boundaryType), boundaryType, null);
        }
        return a;
    }

    public double[,,] CalculateProbabilityDensity()
    {
        var probabilityDensity = new double[Dimensions.x, Dimensions.y, Dimensions.z];

        for (var x = 0; x < Dimensions.x; x++)
        {
            for (var y = 0; y < Dimensions.y; y++)
            {
                for (var z = 0; z < Dimensions.z; z++)
                {
                    var psi = _wavefunction[x, y, z];
                    probabilityDensity[x, y, z] = psi.Magnitude * psi.Magnitude; // |psi|^2
                }
            }
        }

        return probabilityDensity;
    }

    public void InitializeGaussianPacket(double x0, double y0, double z0, double sigma, double kx, double ky, double kz)
    {
        var normalizedAmplitude = CalculateNormalizationConstant(sigma);

        for (var x = 0; x < Dimensions.x; x++)
        {
            for (var y = 0; y < Dimensions.y; y++)
            {
                for (var z = 0; z < Dimensions.z; z++)
                {
                    var exponent = -((x - x0) * (x - x0) + (y - y0) * (y - y0) + (z - z0) * (z - z0)) / (4 * sigma * sigma);
                    var phase = kx * x + ky * y + kz * z;

                    var realPart = Math.Exp(exponent) * Math.Cos(phase);
                    var imaginaryPart = Math.Exp(exponent) * Math.Sin(phase);

                    _wavefunction[x, y, z] = new Complex(realPart, imaginaryPart) * normalizedAmplitude;
                    _laplacianWavefunction[x, y, z] = CalculateLaplacian(_wavefunction, x, y, z, _boundaryType);
                }
            }
        }
    }

    private double CalculateNormalizationConstant(double sigma)
    {
        var sum = 0.0;

        for (var x = 0; x < Dimensions.x; x++)
        {
            for (var y = 0; y < Dimensions.y; y++)
            {
                for (var z = 0; z < Dimensions.z; z++)
                {
                    var temp = Math.Exp(-((x - Dimensions.x / 2.0) * (x - Dimensions.x / 2.0) +
                                          (y - Dimensions.y / 2.0) * (y - Dimensions.y / 2.0) +
                                          (z - Dimensions.z / 2.0) * (z - Dimensions.z / 2.0)) / (4 * sigma * sigma));
                    sum += temp * temp;
                }
            }
        }

        return 1.0 / Math.Sqrt(sum);
    }

    public float[] GetProbabilityData()
    {
        var probabilityDensity = CalculateProbabilityDensity();
        var flattenedData = new List<float>();

        for (var x = 0; x < Dimensions.x; x++)
        {
            for (var y = 0; y < Dimensions.y; y++)
            {
                for (var z = 0; z < Dimensions.z; z++)
                {
                    var probability = (float)probabilityDensity[x, y, z];

                    // Check for invalid numbers
                    if (float.IsNaN(probability) || float.IsInfinity(probability))
                    {
                        //Logger.LogWarning($"Invalid probability value detected at ({x}, {y}, {z}): {probability}");
                        probability = 0;  // Or handle this case as appropriate
                    }

                    flattenedData.Add(probability);
                }
            }
        }

        // Optional: Normalize the data
        var maxProbability = flattenedData.Any() ? flattenedData.Max() : 0.0f;
        for (var i = 0; i < flattenedData.Count; i++)
        {
            flattenedData[i] /= maxProbability;
        }

        return flattenedData.ToArray();
    }
}
//public Complex GetWavefunctionValue(int x, int y, int z)
//    {
//        if (x < 0 || x >= Dimensions.x || y < 0 || y >= Dimensions.y || z < 0 || z >= Dimensions.z)
//        {
//            throw new ArgumentException("Invalid coordinates.");
//        }

//        return wavefunction[x, y, z];
//    }
//public void InitializeWavefunction(Complex initialWavefunction)
//{
//    for (int x = 0; x < dimensions.x; x++)
//    {
//        for (int y = 0; y < dimensions.y; y++)
//        {
//            for (int z = 0; z < dimensions.z; z++)
//            {
//                wavefunction[x, y, z] = initialWavefunction;
//            }
//        }
//    }
//}

//public void InitializeParticleInBoxWavefunction(int n_x, int n_y, int n_z)
//{
//    wavefunction = GetParticleInABoxWavefunction(n_x, n_y, n_z);
//}

//private Complex[,,] GetParticleInABoxWavefunction(int nX, int nY, int nZ)
//{
//    if (dimensions.x <= 0 || dimensions.y <= 0 || dimensions.z <= 0)
//    {
//        throw new ArgumentException("Invalid box dimensions");
//    }

//    if (nX < 1 || nY < 1 || nZ < 1)
//    {
//        throw new ArgumentException("Invalid quantum numbers");
//    }

//    double Lx = dimensions.x * deltaX;
//    double Ly = dimensions.y * deltaY;
//    double Lz = dimensions.z * deltaZ;

//    var boxWavefunction = new Complex[dimensions.x, dimensions.y, dimensions.z];

//    for (int x = 0; x < dimensions.x; x++)
//    {
//        for (int y = 0; y < dimensions.y; y++)
//        {
//            for (int z = 0; z < dimensions.z; z++)
//            {
//                double psi_x = Math.Sqrt(2.0 / Lx) * Math.Sin((nX * Math.PI * (x * deltaX)) / Lx);
//                double psi_y = Math.Sqrt(2.0 / Ly) * Math.Sin((nY * Math.PI * (y * deltaY)) / Ly);
//                double psi_z = Math.Sqrt(2.0 / Lz) * Math.Sin((nZ * Math.PI * (z * deltaZ)) / Lz);

//                boxWavefunction[x, y, z] = new Complex(psi_x * psi_y * psi_z, 0.0);
//            }
//        }
//    }

//    return boxWavefunction;
//}
//public Complex[,,] GetWavefunction()
//{
//    return wavefunction;
//}

//public double MeasureObservable(Complex[,,] observableOperator)
//{
//    double observableValue = 0.0;

//    // Implement the calculation for the expectation value of the observable
//    // This typically involves integrating the product of the wavefunction,
//    // the observable operator, and the complex conjugate of the wavefunction

//    return observableValue;
//}
//public void UpdatePotential(double[,,] newPotential)
//{
//    // Basic validation (optional)
//    if (newPotential.GetLength(0) != Dimensions.x ||
//        newPotential.GetLength(1) != Dimensions.y ||
//        newPotential.GetLength(2) != Dimensions.z)
//    {
//        throw new ArgumentException("New potential dimensions do not match the system dimensions.");
//    }

//    this.potential = newPotential;
//}

//public void ApplySpecificPotentialChange(PotentialChangeType type, double parameter)
//{
//    switch (type)
//    {
//        case PotentialChangeType.Well:
//            ApplyWellPotential(parameter);
//            break;
//        case PotentialChangeType.Barrier:
//            ApplyBarrierPotential(parameter);
//            break;
//        case PotentialChangeType.ExternalField:
//            ApplyExternalField(parameter);
//            break;
//        case PotentialChangeType.RandomPerturbation:
//            ApplyRandomPerturbation(parameter);
//            break;
//    }
//}

//private void ApplyWellPotential(double depth)
//{
//    // Implementation to create a potential well with the specified depth
//}

//private void ApplyBarrierPotential(double height)
//{
//    // Implementation to create a potential barrier with the specified height
//}
//private void ApplyExternalField(double fieldStrength)
//{
//    // Example: Apply a linear potential that simulates an external field
//    for (int x = 0; x < Dimensions.x; x++)
//    {
//        for (int y = 0; y < Dimensions.y; y++)
//        {
//            for (int z = 0; z < Dimensions.z; z++)
//            {
//                potential[x, y, z] += fieldStrength * x; // Linear potential
//            }
//        }
//    }
//}

//private void ApplyRandomPerturbation(double maxPerturbation)
//{
//    Random random = new Random();
//    for (int x = 0; x < Dimensions.x; x++)
//    {
//        for (int y = 0; y < Dimensions.y; y++)
//        {
//            for (int z = 0; z < Dimensions.z; z++)
//            {
//                // Add a random value to the potential at each point
//                potential[x, y, z] += (random.NextDouble() * 2 - 1) * maxPerturbation;
//            }
//        }
//    }
//}

//private void InitializeSinusoidal(double[] parameters)
//{
//    // Extract parameters like nx, ny, nz, etc.
//    // Fill wavefunction array with sinusoidal wavefunction values
//}