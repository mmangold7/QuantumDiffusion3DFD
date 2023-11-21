using MathNet.Numerics.IntegralTransforms;
using System.Diagnostics;
using System.Numerics;

namespace QuantumDiffusion3DFD.Shared;

public class QuantumSystem
{
    // todo
    // 

    private readonly int _xDimension;
    private readonly int _yDimension;
    private readonly int _zDimension;
    private Complex[]? _sliceX, _sliceY, _sliceZ; // Initialize these arrays once with appropriate sizes
    private readonly BoundaryType _boundaryType;
    private Complex[,,] _wavefunction;
    private readonly Complex[,,] _laplacianWavefunction;
    private readonly double[,,] _potential;
    private readonly double _hbar;
    private readonly double _mass;
    private readonly double _timeStep;
    private readonly double _deltaX;

    public QuantumSystem(
        int x,
        int y,
        int z,
        BoundaryType boundaryType,
        double timeStep, 
        double spaceStep, 
        double mass,
        double hbar)
    {
        _xDimension = x;
        _yDimension = y;
        _zDimension = z;

        _boundaryType = boundaryType;
        _deltaX = spaceStep;
        _timeStep = timeStep;

        _wavefunction = new Complex[x, y, z];
        _laplacianWavefunction = new Complex[x, y, z];
        _potential = new double[x, y, z];

        _mass = mass;
        _hbar = hbar;

        _sliceX = new Complex[x];
        _sliceY = new Complex[y];
        _sliceZ = new Complex[z];
    }

    public void InitializeGaussianPacket(double x0, double y0, double z0, double sigma, double kx, double ky, double kz)
    {
        var normalizedAmplitude = CalculateNormalizationConstant(sigma);

        for (var x = 0; x < _xDimension; x++)
        {
            for (var y = 0; y < _yDimension; y++)
            {
                for (var z = 0; z < _zDimension; z++)
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

    public void ApplySingleTimeEvolutionStepEuler()
    {
        var newWavefunction = new Complex[_xDimension, _yDimension, _zDimension];

        for (var x = 0; x < _xDimension; x++)
        {
            for (var y = 0; y < _yDimension; y++)
            {
                for (var z = 0; z < _zDimension; z++)
                {
                    _laplacianWavefunction[x, y, z] = CalculateLaplacian(_wavefunction, x, y, z, _boundaryType);

                    var timeDerivative = (-_hbar * _hbar / (2 * _mass)) * _laplacianWavefunction[x, y, z] + _potential[x, y, z] * _wavefunction[x, y, z];
                    newWavefunction[x, y, z] = _wavefunction[x, y, z] - (Complex.ImaginaryOne / _hbar) * timeDerivative * _timeStep;
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
                plusIndex = Mod(x + 1, _xDimension, boundaryType);
                minusIndex = Mod(x - 1, _xDimension, boundaryType);
                valuePlus = plusIndex != -1 ? wavefunction[plusIndex, y, z] : 0;
                valueMinus = minusIndex != -1 ? wavefunction[minusIndex, y, z] : 0;
                break;
            case 1: // y-dimension
                plusIndex = Mod(y + 1, _yDimension, boundaryType);
                minusIndex = Mod(y - 1, _yDimension, boundaryType);
                valuePlus = plusIndex != -1 ? wavefunction[x, plusIndex, z] : 0;
                valueMinus = minusIndex != -1 ? wavefunction[x, minusIndex, z] : 0;
                break;
            case 2: // z-dimension
                plusIndex = Mod(z + 1, _zDimension, boundaryType);
                minusIndex = Mod(z - 1, _zDimension, boundaryType);
                valuePlus = plusIndex != -1 ? wavefunction[x, y, plusIndex] : 0;
                valueMinus = minusIndex != -1 ? wavefunction[x, y, minusIndex] : 0;
                break;
            default:
                throw new ArgumentException("Invalid dimension");
        }

        return (valuePlus - 2 * wavefunction[x, y, z] + valueMinus) / dx2;
    }

    private static int Mod(int index, int max, BoundaryType boundaryType)
    {
        switch (boundaryType)
        {
            case BoundaryType.Reflective:
                // For reflective boundaries, return -1 if the index is outside the range.
                // This signals the calling function to treat these points as having a wavefunction value of zero.
                if (index < 0 || index >= max)
                    return -1;
                return index;

            case BoundaryType.Periodic:
                // For periodic boundaries, wrap around the index.
                if (index < 0)
                    return max + index;
                if (index >= max)
                    return index % max;
                return index;

            case BoundaryType.Absorbing:
                // For absorbing boundaries, similar to reflective but might have different handling.
                if (index < 0 || index >= max)
                    return -1;
                return index;

            default:
                throw new ArgumentOutOfRangeException(nameof(boundaryType), boundaryType, "Invalid boundary type.");
        }
    }

    public float CalculateTotalEnergy()
    {
        var totalEnergy = 0.0;

        for (var x = 0; x < _xDimension; x++)
        {
            for (var y = 0; y < _yDimension; y++)
            {
                for (var z = 0; z < _zDimension; z++)
                {
                    var psi = _wavefunction[x, y, z];
                    var probabilityDensity = psi.Magnitude * psi.Magnitude;

                    var laplacianPsi = _laplacianWavefunction[x, y, z]; // Use the stored Laplacian
                    var kineticEnergy = -(_hbar * _hbar / (2 * _mass)) * (laplacianPsi * Complex.Conjugate(psi)).Real;

                    var potentialEnergy = _potential[x, y, z] * probabilityDensity;
                    totalEnergy += kineticEnergy + potentialEnergy;
                }
            }
        }

        return (float)totalEnergy;
    }

    public double[,,] CalculateProbabilityDensity()
    {
        var probabilityDensity = new double[_xDimension, _yDimension, _zDimension];

        for (var x = 0; x < _xDimension; x++)
        {
            for (var y = 0; y < _yDimension; y++)
            {
                for (var z = 0; z < _zDimension; z++)
                {
                    var psi = _wavefunction[x, y, z];
                    probabilityDensity[x, y, z] = psi.Magnitude * psi.Magnitude; // |psi|^2
                }
            }
        }

        return probabilityDensity;
    }

    private double CalculateNormalizationConstant(double sigma)
    {
        var sum = 0.0;

        for (var x = 0; x < _xDimension; x++)
        {
            for (var y = 0; y < _yDimension; y++)
            {
                for (var z = 0; z < _zDimension; z++)
                {
                    var temp = Math.Exp(-((x - _xDimension / 2.0) * (x - _xDimension / 2.0) +
                                          (y - _yDimension / 2.0) * (y - _yDimension / 2.0) +
                                          (z - _zDimension / 2.0) * (z - _zDimension / 2.0)) / (4 * sigma * sigma));
                    sum += temp * temp;
                }
            }
        }

        return 1.0 / Math.Sqrt(sum);
    }

    public float[] GetNormalizedProbabilityData(bool normalize = true)
    {
        var probabilityDensity = CalculateProbabilityDensity();
        var flattenedData = new List<float>();

        for (var x = 0; x < _xDimension; x++)
        {
            for (var y = 0; y < _yDimension; y++)
            {
                for (var z = 0; z < _zDimension; z++)
                {
                    var probability = (float)probabilityDensity[x, y, z];

                    if (float.IsNaN(probability) || float.IsInfinity(probability))
                        probability = 0;  // Or handle this case as appropriate

                    flattenedData.Add(probability);
                }
            }
        }

        if (normalize)
        {
            var maxProbability = flattenedData.Any() ? flattenedData.Max() : 0.0f;
            for (var i = 0; i < flattenedData.Count; i++)
                flattenedData[i] /= maxProbability;
        }

        return flattenedData.ToArray();
    }

    public void ApplySingleTimeEvolutionStepRK4()
    {
        var newWavefunction = new Complex[_xDimension, _yDimension, _zDimension];

        for (var x = 0; x < _xDimension; x++)
        {
            for (var y = 0; y < _yDimension; y++)
            {
                for (var z = 0; z < _zDimension; z++)
                {
                    // Calculate k1
                    var k1 = TimeDerivative(_wavefunction, x, y, z);

                    // Create intermediate wavefunctions for k2, k3, k4 calculations
                    var wavefunctionK2 = new Complex[_xDimension, _yDimension, _zDimension];
                    var wavefunctionK3 = new Complex[_xDimension, _yDimension, _zDimension];
                    var wavefunctionK4 = new Complex[_xDimension, _yDimension, _zDimension];

                    // Copy current wavefunction and apply k1 for k2 calculation
                    Array.Copy(_wavefunction, wavefunctionK2, _wavefunction.Length);
                    wavefunctionK2[x, y, z] += 0.5 * _timeStep * k1;

                    // Calculate k2
                    var k2 = TimeDerivative(wavefunctionK2, x, y, z);

                    // Copy current wavefunction and apply k2 for k3 calculation
                    Array.Copy(_wavefunction, wavefunctionK3, _wavefunction.Length);
                    wavefunctionK3[x, y, z] += 0.5 * _timeStep * k2;

                    // Calculate k3
                    var k3 = TimeDerivative(wavefunctionK3, x, y, z);

                    // Copy current wavefunction and apply k3 for k4 calculation
                    Array.Copy(_wavefunction, wavefunctionK4, _wavefunction.Length);
                    wavefunctionK4[x, y, z] += _timeStep * k3;

                    // Calculate k4
                    var k4 = TimeDerivative(wavefunctionK4, x, y, z);

                    // Combine k values to update the wavefunction
                    newWavefunction[x, y, z] = _wavefunction[x, y, z] + (_timeStep / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
                }
            }
        }

        _wavefunction = newWavefunction;
    }

    public void ApplySingleTimeEvolutionStepSSFM()
    {
        Profiling.RunWithClockingLog(() => FourierTransform(ref _wavefunction));
        Profiling.RunWithClockingLog(() => KineticEvolution(ref _wavefunction, _timeStep / 2));
        Profiling.RunWithClockingLog(() => InverseFourierTransform(ref _wavefunction));
        Profiling.RunWithClockingLog(() => PotentialEvolution(ref _wavefunction, _timeStep));
        Profiling.RunWithClockingLog(() => FourierTransform(ref _wavefunction));
        Profiling.RunWithClockingLog(() => KineticEvolution(ref _wavefunction, _timeStep / 2));
        Profiling.RunWithClockingLog(() => InverseFourierTransform(ref _wavefunction));
    }

    private Complex TimeDerivative(Complex[,,] wavefunction, int x, int y, int z)
    {
        var laplacian = CalculateLaplacian(wavefunction, x, y, z, _boundaryType);
        return (-_hbar * _hbar / (2 * _mass)) * laplacian + _potential[x, y, z] * wavefunction[x, y, z];
    }

    private void FourierTransform(ref Complex[,,] wavefunction) => Transform3D(wavefunction, Fourier.Forward);

    private void InverseFourierTransform(ref Complex[,,] wavefunction) => Transform3D(wavefunction, Fourier.Inverse);

    private void Transform3D(Complex[,,] data, Action<Complex[]> transform)
    {
        // Transform along x
        for (int y = 0; y < _yDimension; y++)
        {
            for (int z = 0; z < _zDimension; z++)
            {
                for (int x = 0; x < _xDimension; x++)
                {
                    _sliceX[x] = data[x, y, z];
                }

                transform(_sliceX);

                for (int x = 0; x < _xDimension; x++)
                {
                    data[x, y, z] = _sliceX[x];
                }
            }
        }

        // Transform along y
        for (int x = 0; x < _xDimension; x++)
        {
            for (int z = 0; z < _zDimension; z++)
            {
                for (int y = 0; y < _yDimension; y++)
                {
                    _sliceY[y] = data[x, y, z];
                }

                transform(_sliceY);

                for (int y = 0; y < _yDimension; y++)
                {
                    data[x, y, z] = _sliceY[y];
                }
            }
        }

        // Transform along z
        for (int x = 0; x < _xDimension; x++)
        {
            for (int y = 0; y < _yDimension; y++)
            {
                for (int z = 0; z < _zDimension; z++)
                {
                    _sliceZ[z] = data[x, y, z];
                }

                transform(_sliceZ);

                for (int z = 0; z < _zDimension; z++)
                {
                    data[x, y, z] = _sliceZ[z];
                }
            }
        }
    }

    private void KineticEvolution(ref Complex[,,] wavefunction, double timeStep)
    {
        double dx = _deltaX; // Space step
        double dkx = 2 * Math.PI / (_xDimension * dx); // Momentum step in x
        double dky = 2 * Math.PI / (_yDimension * dx); // Assuming same space step in y
        double dkz = 2 * Math.PI / (_zDimension * dx); // Assuming same space step in z

        // Precompute phase shifts if feasible
        // For example, if you know the range of kineticEnergy values
        Dictionary<double, (double cos, double sin)> precomputedValues = new Dictionary<double, (double cos, double sin)>();

        for (int x = 0; x < _xDimension; x++)
        {
            for (int y = 0; y < _yDimension; y++)
            {
                for (int z = 0; z < _zDimension; z++)
                {
                    // Calculate the momentum for each dimension
                    double kx = (x <= _xDimension / 2) ? x * dkx : (x - _xDimension) * dkx;
                    double ky = (y <= _yDimension / 2) ? y * dky : (y - _yDimension) * dky;
                    double kz = (z <= _zDimension / 2) ? z * dkz : (z - _zDimension) * dkz;

                    // Calculate the kinetic energy term
                    double kineticEnergy = (kx * kx + ky * ky + kz * kz) * (_hbar * _hbar) / (2 * _mass);

                    // Update the wavefunction based on the kinetic energy
                    double phaseShift = -kineticEnergy * timeStep / _hbar;

                    if (!precomputedValues.TryGetValue(phaseShift, out var trigValues))
                    {
                        trigValues = (Math.Cos(phaseShift), Math.Sin(phaseShift));
                        precomputedValues[phaseShift] = trigValues;
                    }

                    wavefunction[x, y, z] *= new Complex(trigValues.cos, trigValues.sin);
                }
            }
        }
    }

    private void PotentialEvolution(ref Complex[,,] wavefunction, double timeStep)
    {
        for (int x = 0; x < _xDimension; x++)
        {
            for (int y = 0; y < _yDimension; y++)
            {
                for (int z = 0; z < _zDimension; z++)
                {
                    // Get the potential energy at this point
                    double potentialEnergy = _potential[x, y, z];

                    // Calculate the phase shift due to the potential energy
                    double phaseShift = -potentialEnergy * timeStep / _hbar;

                    // Update the wavefunction at this point
                    wavefunction[x, y, z] *= new Complex(Math.Cos(phaseShift), Math.Sin(phaseShift));
                }
            }
        }
    }
}

//string dll1Path = Path.Combine(System.Environment.CurrentDirectory, "_framework", "libfftw3-3-x64.dll");
//string dll2Path = Path.Combine(System.Environment.CurrentDirectory, "_framework", "libfftw3-3-x86.dll");
//string dll3Path = Path.Combine(System.Environment.CurrentDirectory, "_framework", "libfftw3-3-32.dll");
//string dll4Path = Path.Combine(System.Environment.CurrentDirectory, "_framework", "libfftw3-3-64.dll");
//NativeLibrary.Load(dll1Path);
//NativeLibrary.Load(dll2Path);
//NativeLibrary.Load(dll3Path);
//NativeLibrary.Load(dll4Path);
//FftwInterop.fftw_init_threads();

//private void FourierTransform(ref Complex[,,] wavefunction)
//{
//    var output = new Complex[_xDimension, _yDimension, _zDimension];
//    using var pinIn = new PinnedArray<Complex>(wavefunction);
//    using var pinOut = new PinnedArray<Complex>(output);
//    using var fftPlan = FftwPlanC2C.Create(pinIn, pinOut, DftDirection.Forwards);
//    fftPlan.Execute();
//    Array.Copy(output, wavefunction, wavefunction.Length);
//}

//private void InverseFourierTransform(ref Complex[,,] wavefunction)
//{
//    var output = new Complex[_xDimension, _yDimension, _zDimension];
//    using var pinIn = new PinnedArray<Complex>(wavefunction);
//    using var pinOut = new PinnedArray<Complex>(output);
//    using var ifftPlan = FftwPlanC2C.Create(pinIn, pinOut, DftDirection.Backwards);
//    ifftPlan.Execute();
//    Array.Copy(output, wavefunction, wavefunction.Length);
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