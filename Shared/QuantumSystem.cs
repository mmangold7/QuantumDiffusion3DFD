using MathNet.Numerics.IntegralTransforms;
using System.Diagnostics;
using System.Numerics;

namespace QuantumDiffusion3DFD.Shared;

public class QuantumSystem
{
    private Complex[]? _sliceX, _sliceY, _sliceZ; // Initialize these arrays once with appropriate sizes
    private readonly (int x, int y, int z) _dimensions;
    private readonly BoundaryType _boundaryType;
    private Complex[,,] _wavefunction;
    private readonly Complex[,,] _laplacianWavefunction;
    private readonly double[,,] _potential;
    private readonly double _hbar;
    private readonly double _mass;
    private readonly double _timeStep;
    private readonly double _deltaX;
    Stopwatch? stopwatch;

    public QuantumSystem(
        (int x, int y, int z) dimensions,
        BoundaryType boundaryType,
        double timeStep, 
        double spaceStep, 
        double mass,
        double hbar)
    {
        _dimensions = dimensions;
        _boundaryType = boundaryType;
        _deltaX = spaceStep;
        _timeStep = timeStep;
        _wavefunction = new Complex[dimensions.x, dimensions.y, dimensions.z];
        _laplacianWavefunction = new Complex[dimensions.x, dimensions.y, dimensions.z];
        _potential = new double[dimensions.x, dimensions.y, dimensions.z];
        _mass = mass;
        _hbar = hbar;

        _sliceX = new Complex[dimensions.x];
        _sliceY = new Complex[dimensions.y];
        _sliceZ = new Complex[dimensions.z];
    }

    public double CalculateTotalEnergy()
    {
        var totalEnergy = 0.0;

        for (var x = 0; x < _dimensions.x; x++)
        {
            for (var y = 0; y < _dimensions.y; y++)
            {
                for (var z = 0; z < _dimensions.z; z++)
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

        return totalEnergy;
    }

    public void ApplySingleTimeEvolutionStepEuler()
    {
        var newWavefunction = new Complex[_dimensions.x, _dimensions.y, _dimensions.z];

        for (var x = 0; x < _dimensions.x; x++)
        {
            for (var y = 0; y < _dimensions.y; y++)
            {
                for (var z = 0; z < _dimensions.z; z++)
                {
                    _laplacianWavefunction[x, y, z] = CalculateLaplacian(_wavefunction, x, y, z, _boundaryType);

                    var timeDerivative = (-_hbar * _hbar / (2 * _mass)) * _laplacianWavefunction[x, y, z] + _potential[x, y, z] * _wavefunction[x, y, z];
                    newWavefunction[x, y, z] = _wavefunction[x, y, z] - (Complex.ImaginaryOne / _hbar) * timeDerivative * _timeStep;
                }
            }
        }
        _wavefunction = newWavefunction;
    }

    public void ApplySingleTimeEvolutionStepRK4()
    {
        var newWavefunction = new Complex[_dimensions.x, _dimensions.y, _dimensions.z];

        for (var x = 0; x < _dimensions.x; x++)
        {
            for (var y = 0; y < _dimensions.y; y++)
            {
                for (var z = 0; z < _dimensions.z; z++)
                {
                    // Calculate k1
                    var k1 = TimeDerivative(_wavefunction, x, y, z);

                    // Create intermediate wavefunctions for k2, k3, k4 calculations
                    var wavefunctionK2 = new Complex[_dimensions.x, _dimensions.y, _dimensions.z];
                    var wavefunctionK3 = new Complex[_dimensions.x, _dimensions.y, _dimensions.z];
                    var wavefunctionK4 = new Complex[_dimensions.x, _dimensions.y, _dimensions.z];

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
        //string dll1Path = Path.Combine(System.Environment.CurrentDirectory, "_framework", "libfftw3-3-x64.dll");
        //string dll2Path = Path.Combine(System.Environment.CurrentDirectory, "_framework", "libfftw3-3-x86.dll");
        //string dll3Path = Path.Combine(System.Environment.CurrentDirectory, "_framework", "libfftw3-3-32.dll");
        //string dll4Path = Path.Combine(System.Environment.CurrentDirectory, "_framework", "libfftw3-3-64.dll");

        //NativeLibrary.Load(dll1Path);
        //NativeLibrary.Load(dll2Path);
        //NativeLibrary.Load(dll3Path);
        //NativeLibrary.Load(dll4Path);

        //FftwInterop.fftw_init_threads();

        var debugOutput = false;
        if (debugOutput) stopwatch = Stopwatch.StartNew();

        FourierTransform(ref _wavefunction);
        if (debugOutput) Extensions.LogMethodTime(nameof(FourierTransform), stopwatch);

        KineticEvolution(ref _wavefunction, _timeStep / 2);
        if (debugOutput) Extensions.LogMethodTime(nameof(KineticEvolution), stopwatch);

        InverseFourierTransform(ref _wavefunction);
        if (debugOutput) Extensions.LogMethodTime(nameof(InverseFourierTransform), stopwatch);

        PotentialEvolution(ref _wavefunction, _timeStep);
        if (debugOutput) Extensions.LogMethodTime(nameof(PotentialEvolution), stopwatch);

        FourierTransform(ref _wavefunction);
        if (debugOutput) Extensions.LogMethodTime(nameof(FourierTransform), stopwatch);

        KineticEvolution(ref _wavefunction, _timeStep / 2);
        if (debugOutput) Extensions.LogMethodTime(nameof(KineticEvolution), stopwatch);

        InverseFourierTransform(ref _wavefunction);
        if (debugOutput) Extensions.LogMethodTime(nameof(InverseFourierTransform), stopwatch);
    }

    private Complex TimeDerivative(Complex[,,] wavefunction, int x, int y, int z)
    {
        var laplacian = CalculateLaplacian(wavefunction, x, y, z, _boundaryType);
        return (-_hbar * _hbar / (2 * _mass)) * laplacian + _potential[x, y, z] * wavefunction[x, y, z];
    }

    //private void FourierTransform(ref Complex[,,] wavefunction)
    //{
    //    var output = new Complex[_dimensions.x, _dimensions.y, _dimensions.z];
    //    using var pinIn = new PinnedArray<Complex>(wavefunction);
    //    using var pinOut = new PinnedArray<Complex>(output);
    //    using var fftPlan = FftwPlanC2C.Create(pinIn, pinOut, DftDirection.Forwards);
    //    fftPlan.Execute();
    //    Array.Copy(output, wavefunction, wavefunction.Length);
    //}

    //private void InverseFourierTransform(ref Complex[,,] wavefunction)
    //{
    //    var output = new Complex[_dimensions.x, _dimensions.y, _dimensions.z];
    //    using var pinIn = new PinnedArray<Complex>(wavefunction);
    //    using var pinOut = new PinnedArray<Complex>(output);
    //    using var ifftPlan = FftwPlanC2C.Create(pinIn, pinOut, DftDirection.Backwards);
    //    ifftPlan.Execute();
    //    Array.Copy(output, wavefunction, wavefunction.Length);
    //}

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
                plusIndex = Mod(x + 1, _dimensions.x, boundaryType);
                minusIndex = Mod(x - 1, _dimensions.x, boundaryType);
                valuePlus = plusIndex != -1 ? wavefunction[plusIndex, y, z] : 0;
                valueMinus = minusIndex != -1 ? wavefunction[minusIndex, y, z] : 0;
                break;
            case 1: // y-dimension
                plusIndex = Mod(y + 1, _dimensions.y, boundaryType);
                minusIndex = Mod(y - 1, _dimensions.y, boundaryType);
                valuePlus = plusIndex != -1 ? wavefunction[x, plusIndex, z] : 0;
                valueMinus = minusIndex != -1 ? wavefunction[x, minusIndex, z] : 0;
                break;
            case 2: // z-dimension
                plusIndex = Mod(z + 1, _dimensions.z, boundaryType);
                minusIndex = Mod(z - 1, _dimensions.z, boundaryType);
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
        var probabilityDensity = new double[_dimensions.x, _dimensions.y, _dimensions.z];

        for (var x = 0; x < _dimensions.x; x++)
        {
            for (var y = 0; y < _dimensions.y; y++)
            {
                for (var z = 0; z < _dimensions.z; z++)
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

        for (var x = 0; x < _dimensions.x; x++)
        {
            for (var y = 0; y < _dimensions.y; y++)
            {
                for (var z = 0; z < _dimensions.z; z++)
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

        for (var x = 0; x < _dimensions.x; x++)
        {
            for (var y = 0; y < _dimensions.y; y++)
            {
                for (var z = 0; z < _dimensions.z; z++)
                {
                    var temp = Math.Exp(-((x - _dimensions.x / 2.0) * (x - _dimensions.x / 2.0) +
                                          (y - _dimensions.y / 2.0) * (y - _dimensions.y / 2.0) +
                                          (z - _dimensions.z / 2.0) * (z - _dimensions.z / 2.0)) / (4 * sigma * sigma));
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

        for (var x = 0; x < _dimensions.x; x++)
        {
            for (var y = 0; y < _dimensions.y; y++)
            {
                for (var z = 0; z < _dimensions.z; z++)
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

    private void FourierTransform(ref Complex[,,] wavefunction) => Transform3D(wavefunction, Fourier.Forward);

    private void InverseFourierTransform(ref Complex[,,] wavefunction) => Transform3D(wavefunction, Fourier.Inverse);

    private void Transform3D(Complex[,,] data, Action<Complex[]> transform)
    {
        // Transform along x
        for (int y = 0; y < _dimensions.y; y++)
        {
            for (int z = 0; z < _dimensions.z; z++)
            {
                for (int x = 0; x < _dimensions.x; x++)
                {
                    _sliceX[x] = data[x, y, z];
                }

                transform(_sliceX);

                for (int x = 0; x < _dimensions.x; x++)
                {
                    data[x, y, z] = _sliceX[x];
                }
            }
        }

        // Transform along y
        for (int x = 0; x < _dimensions.x; x++)
        {
            for (int z = 0; z < _dimensions.z; z++)
            {
                for (int y = 0; y < _dimensions.y; y++)
                {
                    _sliceY[y] = data[x, y, z];
                }

                transform(_sliceY);

                for (int y = 0; y < _dimensions.y; y++)
                {
                    data[x, y, z] = _sliceY[y];
                }
            }
        }

        // Transform along z
        for (int x = 0; x < _dimensions.x; x++)
        {
            for (int y = 0; y < _dimensions.y; y++)
            {
                for (int z = 0; z < _dimensions.z; z++)
                {
                    _sliceZ[z] = data[x, y, z];
                }

                transform(_sliceZ);

                for (int z = 0; z < _dimensions.z; z++)
                {
                    data[x, y, z] = _sliceZ[z];
                }
            }
        }
    }

    private void KineticEvolution(ref Complex[,,] wavefunction, double timeStep)
    {
        double dx = _deltaX; // Space step
        double dkx = 2 * Math.PI / (_dimensions.x * dx); // Momentum step in x
        double dky = 2 * Math.PI / (_dimensions.y * dx); // Assuming same space step in y
        double dkz = 2 * Math.PI / (_dimensions.z * dx); // Assuming same space step in z

        // Precompute phase shifts if feasible
        // For example, if you know the range of kineticEnergy values
        Dictionary<double, (double cos, double sin)> precomputedValues = new Dictionary<double, (double cos, double sin)>();

        for (int x = 0; x < _dimensions.x; x++)
        {
            for (int y = 0; y < _dimensions.y; y++)
            {
                for (int z = 0; z < _dimensions.z; z++)
                {
                    // Calculate the momentum for each dimension
                    double kx = (x <= _dimensions.x / 2) ? x * dkx : (x - _dimensions.x) * dkx;
                    double ky = (y <= _dimensions.y / 2) ? y * dky : (y - _dimensions.y) * dky;
                    double kz = (z <= _dimensions.z / 2) ? z * dkz : (z - _dimensions.z) * dkz;

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
        for (int x = 0; x < _dimensions.x; x++)
        {
            for (int y = 0; y < _dimensions.y; y++)
            {
                for (int z = 0; z < _dimensions.z; z++)
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
//    // Extract parameters like _dimensions.x, _dimensions.y, _dimensions.z, etc.
//    // Fill wavefunction array with sinusoidal wavefunction values
//}