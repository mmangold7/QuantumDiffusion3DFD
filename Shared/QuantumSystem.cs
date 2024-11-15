using System.Numerics;
using QuantumDiffusion3DFD.Shared.Extensions;
using QuantumDiffusion3DFD.Shared.Models;
using QuantumDiffusion3DFD.Shared.Models.Enums;

namespace QuantumDiffusion3DFD.Shared;

public class QuantumSystem
{
    private readonly float _hbar;
    private readonly float _mass;
    private readonly float _timeStep;
    private readonly float _spaceStep;
    private readonly int _xDimension, _yDimension, _zDimension;
    private readonly BoundaryType _boundaryType;

    private Complex[]? _sliceX, _sliceY, _sliceZ;
    private Complex[,,] _wavefunction;
    private readonly Complex[,,] _laplacian;
    private readonly float[,,] _potential;
    private Complex[,,] _absorbingPotential;
    private readonly float[] _probability;
    private readonly float[] _previousProbabilityData;
    private float? _originalTotalEnergy;

    public QuantumSystem(
        int x, int y, int z,
        BoundaryType boundaryType,
        float timeStep, float spaceStep,
        float mass, float hbar)
    {
        _mass = mass;
        _hbar = hbar;
        _spaceStep = spaceStep;
        _timeStep = timeStep;
        _xDimension = x;
        _yDimension = y;
        _zDimension = z;
        _boundaryType = boundaryType;
        _sliceX = new Complex[x];
        _sliceY = new Complex[y];
        _sliceZ = new Complex[z];
        _wavefunction = new Complex[x, y, z];
        _laplacian = new Complex[x, y, z];
        _potential = new float[x, y, z];
        _probability = new float[x * y * z];
        _previousProbabilityData = new float[x * y * z];

        if (boundaryType == BoundaryType.Absorbing)
        {
            InitializeAbsorbingPotential();
        }
    }

    private void InitializeAbsorbingPotential()
    {
        _absorbingPotential = new Complex[_xDimension, _yDimension, _zDimension];
        // Define the functional form of the potential here.
        // For example, a simple polynomial that increases towards the edges.
        for (int x = 0; x < _xDimension; x++)
        {
            for (int y = 0; y < _yDimension; y++)
            {
                for (int z = 0; z < _zDimension; z++)
                {
                    // Example: a polynomial potential that activates near the boundaries
                    _absorbingPotential[x, y, z] = CalculateAbsorbingPotential(x, y, z);
                }
            }
        }
    }

    private Complex CalculateAbsorbingPotential(int x, int y, int z)
    {
        // Get the distance to the nearest boundary
        double distanceToBoundary = CalculateDistanceToBoundary(x, y, z);

        // Implement a stronger absorbing potential based on this distance
        // Example: using a higher power and a scaling factor
        double scaleFactor = 3.0; // Adjust this value as needed
        double potentialStrength = -Math.Pow(distanceToBoundary, 4) * scaleFactor;

        return new Complex(0, potentialStrength);
    }

    private double CalculateDistanceToBoundary(int x, int y, int z)
    {
        // Calculate the distances to the lower and upper boundaries along each axis
        double distanceToXBoundary = Math.Min(x, _xDimension - x);
        double distanceToYBoundary = Math.Min(y, _yDimension - y);
        double distanceToZBoundary = Math.Min(z, _zDimension - z);

        // The distance to the nearest boundary is the minimum of these three distances
        return Math.Min(Math.Min(distanceToXBoundary, distanceToYBoundary), distanceToZBoundary);
    }

    public void InitializeGaussianPacket(
        double x0, double y0, double z0, double sigma, double kx, double ky, double kz)
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
                    _laplacian[x, y, z] = CalculateLaplacian(_wavefunction, x, y, z, _boundaryType);
                }
            }
        }
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

    public QuantumState UpdateSimulation(bool onlySignificantlyChanged)
    {
        var newTotalEnergy = 0.0f;
        var flattenedData = new float[_xDimension * _yDimension * _zDimension];

        for (var x = 0; x < _xDimension; x++)
        {
            for (var y = 0; y < _yDimension; y++)
            {
                for (var z = 0; z < _zDimension; z++)
                {
                    var pointWaveValue = UpdateWavefunctionPoint(x, y, z);
                    var psiMag = (float)pointWaveValue.Magnitude;
                    var potentialAtPoint = _potential[x, y, z];
                    var laplacianAtPoint = _laplacian[x, y, z];
                    var pointProbabilityDensity = psiMag * psiMag;
                    var flatIndex = x * _yDimension * _zDimension + y * _zDimension + z;
                    _probability[flatIndex] = pointProbabilityDensity;
                    newTotalEnergy += CalculateEnergyAtPoint(
                        pointWaveValue, pointProbabilityDensity, potentialAtPoint, laplacianAtPoint);
                }
            }
        }

        var newMaxProbability = _probability.Max();
        for (var i = 0; i < _probability.Length; i++)
            flattenedData[i] = _probability[i] / newMaxProbability;

        _originalTotalEnergy ??= newTotalEnergy;

        return new QuantumState
        {
            OriginalTotalEnergy = _originalTotalEnergy.Value,
            CurrentTotalEnergy = newTotalEnergy,
            ProbabilityData = GetProbability(flattenedData, onlySignificantlyChanged)
        };
    }

    private Complex UpdateWavefunctionPoint(int x, int y, int z)
    {
        _laplacian[x, y, z] = CalculateLaplacian(_wavefunction, x, y, z, _boundaryType);
        var timeDerivative = (-_hbar * _hbar / (2 * _mass)) * _laplacian[x, y, z]
                             + _potential[x, y, z] * _wavefunction[x, y, z];

        // Apply the absorbing potential only for absorbing boundary type
        if (_boundaryType == BoundaryType.Absorbing)
        {
            timeDerivative += _absorbingPotential[x, y, z] * _wavefunction[x, y, z];
        }

        var newPointValue = _wavefunction[x, y, z] - (Complex.ImaginaryOne / _hbar) * timeDerivative * _timeStep;
        _wavefunction[x, y, z] = newPointValue;
        return newPointValue;
    }

    private float CalculateEnergyAtPoint(
        Complex psi, float probabilityDensityAtPoint, float potentialAtPoint, Complex laplacianAtPoint)
    {
        var potentialEnergy = potentialAtPoint * probabilityDensityAtPoint;
        var kineticEnergy = -(_hbar * _hbar / (2 * _mass)) * (float)((laplacianAtPoint * Complex.Conjugate(psi)).Real);
        var totalEnergy = kineticEnergy + potentialEnergy;
        return totalEnergy;
    }

    public List<object> GetProbability(float[] probabilityData, bool onlySignificantlyChanged = false)
    {
        var updatedData = new List<object>();
        var maxProbability = 1.0f;
        var updateThreshold = maxProbability * 0.1f;
        var opacityScale = 0.75f;

        for (int i = 0; i < probabilityData.Length; i++)
        {
            var newProbability = probabilityData[i];
            if (!onlySignificantlyChanged || Math.Abs(newProbability - _previousProbabilityData[i]) > updateThreshold)
            {
                var color = GraphicsExtensions.InterpolateColor(newProbability);
                var opacity = GraphicsExtensions.SigmoidOpacity(newProbability * opacityScale);
                updatedData.Add(new { index = i, color, opacity });
                _previousProbabilityData[i] = newProbability;
            }
        }

        return updatedData;
    }

    private Complex CalculateLaplacian(Complex[,,] wavefunction, int x, int y, int z, BoundaryType boundaryType)
    {
        var dx2 = _spaceStep * _spaceStep;

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
                break;
            case 1: // y-dimension
                plusIndex = Mod(y + 1, _yDimension, boundaryType);
                minusIndex = Mod(y - 1, _yDimension, boundaryType);
                break;
            case 2: // z-dimension
                plusIndex = Mod(z + 1, _zDimension, boundaryType);
                minusIndex = Mod(z - 1, _zDimension, boundaryType);
                break;
            default:
                throw new ArgumentException("Invalid dimension");
        }

        // Handle the boundary conditions
        if (boundaryType == BoundaryType.Absorbing)
        {
            valuePlus = wavefunction[plusIndex, y, z];
            valueMinus = wavefunction[minusIndex, y, z];
        }
        else
        {
            valuePlus = plusIndex != -1 ? wavefunction[plusIndex, y, z] : 0;
            valueMinus = minusIndex != -1 ? wavefunction[minusIndex, y, z] : 0;
        }

        return (valuePlus - 2 * wavefunction[x, y, z] + valueMinus) / dx2;
    }

    private static int Mod(int index, int max, BoundaryType boundaryType)
    {
        switch (boundaryType)
        {
            case BoundaryType.Reflective:
                if (index < 0 || index >= max)
                    return -1;
                return index;

            case BoundaryType.Periodic:
                if (index < 0)
                    return max + index;
                if (index >= max)
                    return index % max;
                return index;

            case BoundaryType.Absorbing:
                return Math.Clamp(index, 0, max - 1); // Ensure index stays within bounds

            default:
                throw new ArgumentOutOfRangeException(nameof(boundaryType), boundaryType, "Invalid boundary type.");
        }
    }
}


//public void ApplySingleTimeEvolutionStepRK4()
//{
//    var newWavefunction = new Complex[_xDimension, _yDimension, _zDimension];

//    for (var x = 0; x < _xDimension; x++)
//    {
//        for (var y = 0; y < _yDimension; y++)
//        {
//            for (var z = 0; z < _zDimension; z++)
//            {
//                // Calculate k1
//                var k1 = TimeDerivative(_wavefunction, x, y, z);

//                // Create intermediate wavefunctions for k2, k3, k4 calculations
//                var wavefunctionK2 = new Complex[_xDimension, _yDimension, _zDimension];
//                var wavefunctionK3 = new Complex[_xDimension, _yDimension, _zDimension];
//                var wavefunctionK4 = new Complex[_xDimension, _yDimension, _zDimension];

//                // Copy current wavefunction and apply k1 for k2 calculation
//                Array.Copy(_wavefunction, wavefunctionK2, _wavefunction.Length);
//                wavefunctionK2[x, y, z] += 0.5 * _timeStep * k1;

//                // Calculate k2
//                var k2 = TimeDerivative(wavefunctionK2, x, y, z);

//                // Copy current wavefunction and apply k2 for k3 calculation
//                Array.Copy(_wavefunction, wavefunctionK3, _wavefunction.Length);
//                wavefunctionK3[x, y, z] += 0.5 * _timeStep * k2;

//                // Calculate k3
//                var k3 = TimeDerivative(wavefunctionK3, x, y, z);

//                // Copy current wavefunction and apply k3 for k4 calculation
//                Array.Copy(_wavefunction, wavefunctionK4, _wavefunction.Length);
//                wavefunctionK4[x, y, z] += _timeStep * k3;

//                // Calculate k4
//                var k4 = TimeDerivative(wavefunctionK4, x, y, z);

//                // Combine k values to update the wavefunction
//                newWavefunction[x, y, z] = _wavefunction[x, y, z] + (_timeStep / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
//            }
//        }
//    }

//    _wavefunction = newWavefunction;
//}

//public void ApplySingleTimeEvolutionStepSSFM()
//{
//    Profiling.RunWithClockingLog(() => FourierTransform(ref _wavefunction));
//    Profiling.RunWithClockingLog(() => KineticEvolution(ref _wavefunction, _timeStep / 2));
//    Profiling.RunWithClockingLog(() => InverseFourierTransform(ref _wavefunction));
//    Profiling.RunWithClockingLog(() => PotentialEvolution(ref _wavefunction, _timeStep));
//    Profiling.RunWithClockingLog(() => FourierTransform(ref _wavefunction));
//    Profiling.RunWithClockingLog(() => KineticEvolution(ref _wavefunction, _timeStep / 2));
//    Profiling.RunWithClockingLog(() => InverseFourierTransform(ref _wavefunction));
//}

//private Complex TimeDerivative(Complex[,,] wavefunction, int x, int y, int z)
//{
//    var laplacian = CalculateLaplacian(wavefunction, x, y, z, _boundaryType);
//    return (-_hbar * _hbar / (2 * _mass)) * laplacian + _potential[x, y, z] * wavefunction[x, y, z];
//}

//private void FourierTransform(ref Complex[,,] wavefunction) => Transform3D(wavefunction, Fourier.Forward);

//private void InverseFourierTransform(ref Complex[,,] wavefunction) => Transform3D(wavefunction, Fourier.Inverse);

//private void Transform3D(Complex[,,] data, Action<Complex[]> transform)
//{
//    // Transform along x
//    for (int y = 0; y < _yDimension; y++)
//    {
//        for (int z = 0; z < _zDimension; z++)
//        {
//            for (int x = 0; x < _xDimension; x++)
//            {
//                _sliceX[x] = data[x, y, z];
//            }

//            transform(_sliceX);

//            for (int x = 0; x < _xDimension; x++)
//            {
//                data[x, y, z] = _sliceX[x];
//            }
//        }
//    }

//    // Transform along y
//    for (int x = 0; x < _xDimension; x++)
//    {
//        for (int z = 0; z < _zDimension; z++)
//        {
//            for (int y = 0; y < _yDimension; y++)
//            {
//                _sliceY[y] = data[x, y, z];
//            }

//            transform(_sliceY);

//            for (int y = 0; y < _yDimension; y++)
//            {
//                data[x, y, z] = _sliceY[y];
//            }
//        }
//    }

//    // Transform along z
//    for (int x = 0; x < _xDimension; x++)
//    {
//        for (int y = 0; y < _yDimension; y++)
//        {
//            for (int z = 0; z < _zDimension; z++)
//            {
//                _sliceZ[z] = data[x, y, z];
//            }

//            transform(_sliceZ);

//            for (int z = 0; z < _zDimension; z++)
//            {
//                data[x, y, z] = _sliceZ[z];
//            }
//        }
//    }
//}

//private void KineticEvolution(ref Complex[,,] wavefunction, double timeStep)
//{
//    double dx = _spaceStep; // Space step
//    double dkx = 2 * Math.PI / (_xDimension * dx); // Momentum step in x
//    double dky = 2 * Math.PI / (_yDimension * dx); // Assuming same space step in y
//    double dkz = 2 * Math.PI / (_zDimension * dx); // Assuming same space step in z

//    // Precompute phase shifts if feasible
//    // For example, if you know the range of kineticEnergy values
//    Dictionary<double, (double cos, double sin)> precomputedValues = new Dictionary<double, (double cos, double sin)>();

//    for (int x = 0; x < _xDimension; x++)
//    {
//        for (int y = 0; y < _yDimension; y++)
//        {
//            for (int z = 0; z < _zDimension; z++)
//            {
//                // Calculate the momentum for each dimension
//                double kx = (x <= _xDimension / 2) ? x * dkx : (x - _xDimension) * dkx;
//                double ky = (y <= _yDimension / 2) ? y * dky : (y - _yDimension) * dky;
//                double kz = (z <= _zDimension / 2) ? z * dkz : (z - _zDimension) * dkz;

//                // Calculate the kinetic energy term
//                double kineticEnergy = (kx * kx + ky * ky + kz * kz) * (_hbar * _hbar) / (2 * _mass);

//                // Update the wavefunction based on the kinetic energy
//                double phaseShift = -kineticEnergy * timeStep / _hbar;

//                if (!precomputedValues.TryGetValue(phaseShift, out var trigValues))
//                {
//                    trigValues = (Math.Cos(phaseShift), Math.Sin(phaseShift));
//                    precomputedValues[phaseShift] = trigValues;
//                }

//                wavefunction[x, y, z] *= new Complex(trigValues.cos, trigValues.sin);
//            }
//        }
//    }
//}

//private void PotentialEvolution(ref Complex[,,] wavefunction, double timeStep)
//{
//    for (int x = 0; x < _xDimension; x++)
//    {
//        for (int y = 0; y < _yDimension; y++)
//        {
//            for (int z = 0; z < _zDimension; z++)
//            {
//                // Get the potential energy at this point
//                double potentialEnergy = _potential[x, y, z];

//                // Calculate the phase shift due to the potential energy
//                double phaseShift = -potentialEnergy * timeStep / _hbar;

//                // Update the wavefunction at this point
//                wavefunction[x, y, z] *= new Complex(Math.Cos(phaseShift), Math.Sin(phaseShift));
//            }
//        }
//    }
//}

//string dll1Path = Path.Combine(System.Environment.CurrentDirectory, "_framework", "libfftw3-3-x64.dll");
//NativeLibrary.Load(dll1Path);

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

//private void ApplyWellPotential(double depth)
//private void ApplyBarrierPotential(double height)
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