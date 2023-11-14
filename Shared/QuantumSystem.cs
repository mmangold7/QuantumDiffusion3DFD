using System.Numerics;

namespace QuantumDiffusion3DFD.Shared;

public class QuantumSystem
{
    private Complex[,,] wavefunction;
    private double[,,] potential;
    private (int x, int y, int z) dimensions;
    private double timeStep;
    private double deltaX, deltaY; // Grid spacing in x and y directions

    // Other necessary constants like particle , Planck's constant, etc.
    private const double hbar = 1.0; // Reduced Planck constant in natural or scaled units
    private const double singleParticleMass = 1.0; // Mass of the particle in scaled/natural units

    private BoundaryType boundaryType; // Field to store the type of boundary condition

    public double TimeStep
    {
        get { return timeStep; }
    }

    public QuantumSystem((int x, int y, int z) dimensions, double timeStep, double spaceStep, BoundaryType boundaryType)
    {
        this.dimensions = dimensions;
        this.timeStep = timeStep;
        this.deltaX = spaceStep; // Set grid spacing in x direction
        this.deltaY = spaceStep; // Set grid spacing in y direction
        this.boundaryType = boundaryType;

        wavefunction = new Complex[dimensions.x, dimensions.y, dimensions.z];
        potential = new double[dimensions.x, dimensions.y, dimensions.z];
        // Initialize wavefunction and potential
    }

    public void InitializeWavefunction(Complex initialWavefunction)
    {
        for (int x = 0; x < dimensions.x; x++)
        {
            for (int y = 0; y < dimensions.y; y++)
            {
                for (int z = 0; z < dimensions.z; z++)
                {
                    wavefunction[x, y, z] = initialWavefunction;
                }
            }
        }
    }

    public void ApplySingleTimeEvolutionStep()
    {
        Complex[,,] newWavefunction = new Complex[dimensions.x, dimensions.y, dimensions.z];

        for (int x = 0; x < dimensions.x; x++)
        {
            for (int y = 0; y < dimensions.y; y++)
            {
                for (int z = 0; z < dimensions.z; z++)
                {
                    // Calculate the Laplacian of psi using the updated method
                    Complex laplacianPsi = CalculateLaplacian(wavefunction, x, y, z, boundaryType);

                    // Schrödinger equation discretized
                    Complex timeDerivative = (-hbar * hbar / (2 * singleParticleMass)) * laplacianPsi + potential[x, y, z] * wavefunction[x, y, z];

                    // Forward Euler method for time stepping
                    newWavefunction[x, y, z] = wavefunction[x, y, z] - (Complex.ImaginaryOne / hbar) * timeDerivative * timeStep;
                }
            }
        }

        wavefunction = newWavefunction;
    }

    private Complex CalculateLaplacian(Complex[,,] wavefunction, int x, int y, int z, BoundaryType boundaryType)
    {
        double dx2 = deltaX * deltaX; // Assuming deltaX is the grid spacing

        Complex d2Psi_dx2 = GetSecondDerivative(wavefunction, x, y, z, 0, boundaryType, dx2);
        Complex d2Psi_dy2 = GetSecondDerivative(wavefunction, x, y, z, 1, boundaryType, dx2);
        Complex d2Psi_dz2 = GetSecondDerivative(wavefunction, x, y, z, 2, boundaryType, dx2);

        return d2Psi_dx2 + d2Psi_dy2 + d2Psi_dz2;
    }

    private Complex GetSecondDerivative(Complex[,,] wavefunction, int x, int y, int z, int dimension, BoundaryType boundaryType, double dx2)
    {
        int plusIndex, minusIndex;
        switch (dimension)
        {
            case 0: // x-dimension
                plusIndex = Mod(x + 1, dimensions.x, boundaryType);
                minusIndex = Mod(x - 1, dimensions.x, boundaryType);
                break;
            case 1: // y-dimension
                plusIndex = Mod(y + 1, dimensions.y, boundaryType);
                minusIndex = Mod(y - 1, dimensions.y, boundaryType);
                break;
            case 2: // z-dimension
                plusIndex = Mod(z + 1, dimensions.z, boundaryType);
                minusIndex = Mod(z - 1, dimensions.z, boundaryType);
                break;
            default:
                throw new ArgumentException("Invalid dimension");
        }

        Complex valuePlus = plusIndex != -1 ? wavefunction[x, y, z] : 0;
        Complex valueMinus = minusIndex != -1 ? wavefunction[x, y, z] : 0;

        return (valuePlus - 2 * wavefunction[x, y, z] + valueMinus) / dx2;
    }

    private int Mod(int a, int b, BoundaryType boundaryType)
    {
        switch (boundaryType)
        {
            case BoundaryType.Reflective:
                if (a < 0 || a >= b)
                    return Math.Abs(b - Math.Abs(a) % b) % b;
                break;
            case BoundaryType.Absorbing:
                if (a < 0 || a >= b)
                    return -1; // Indicates an absorbing boundary
                break;
        }
        return a;
    }

    public void UpdatePotential(double[,,] newPotential)
    {
        // Basic validation (optional)
        if (newPotential.GetLength(0) != dimensions.x ||
            newPotential.GetLength(1) != dimensions.y ||
            newPotential.GetLength(2) != dimensions.z)
        {
            throw new ArgumentException("New potential dimensions do not match the system dimensions.");
        }

        this.potential = newPotential;
    }

    public void ApplySpecificPotentialChange(PotentialChangeType type, double parameter)
    {
        switch (type)
        {
            case PotentialChangeType.Well:
                ApplyWellPotential(parameter);
                break;
            case PotentialChangeType.Barrier:
                ApplyBarrierPotential(parameter);
                break;
            case PotentialChangeType.ExternalField:
                ApplyExternalField(parameter);
                break;
            case PotentialChangeType.RandomPerturbation:
                ApplyRandomPerturbation(parameter);
                break;
        }
    }

    private void ApplyWellPotential(double depth)
    {
        // Implementation to create a potential well with the specified depth
    }

    private void ApplyBarrierPotential(double height)
    {
        // Implementation to create a potential barrier with the specified height
    }

    public double[,,] CalculateProbabilityDensity()
    {
        double[,,] probabilityDensity = new double[dimensions.x, dimensions.y, dimensions.z];

        for (int x = 0; x < dimensions.x; x++)
        {
            for (int y = 0; y < dimensions.y; y++)
            {
                for (int z = 0; z < dimensions.z; z++)
                {
                    Complex psi = wavefunction[x, y, z];
                    probabilityDensity[x, y, z] = psi.Magnitude * psi.Magnitude; // |psi|^2
                }
            }
        }

        return probabilityDensity;
    }

    private void ApplyExternalField(double fieldStrength)
    {
        // Example: Apply a linear potential that simulates an external field
        for (int x = 0; x < dimensions.x; x++)
        {
            for (int y = 0; y < dimensions.y; y++)
            {
                for (int z = 0; z < dimensions.z; z++)
                {
                    potential[x, y, z] += fieldStrength * x; // Linear potential
                }
            }
        }
    }

    private void ApplyRandomPerturbation(double maxPerturbation)
    {
        Random random = new Random();
        for (int x = 0; x < dimensions.x; x++)
        {
            for (int y = 0; y < dimensions.y; y++)
            {
                for (int z = 0; z < dimensions.z; z++)
                {
                    // Add a random value to the potential at each point
                    potential[x, y, z] += (random.NextDouble() * 2 - 1) * maxPerturbation;
                }
            }
        }
    }

    public Complex[,,] GetWavefunction()
    {
        return wavefunction;
    }

    public double MeasureObservable(Complex[,,] observableOperator)
    {
        double observableValue = 0.0;

        // Implement the calculation for the expectation value of the observable
        // This typically involves integrating the product of the wavefunction,
        // the observable operator, and the complex conjugate of the wavefunction

        return observableValue;
    }

    public Complex GetWavefunctionValue(int x, int y, int z)
    {
        if (x < 0 || x >= dimensions.x || y < 0 || y >= dimensions.y || z < 0 || z >= dimensions.z)
        {
            throw new ArgumentException("Invalid coordinates.");
        }

        return wavefunction[x, y, z];
    }
}