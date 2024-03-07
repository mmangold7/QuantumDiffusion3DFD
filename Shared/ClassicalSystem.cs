using System.Numerics;
using QuantumDiffusion3DFD.Shared.Models.Enums;

namespace QuantumDiffusion3DFD.Shared;

public class ClassicalSystem
{
    private readonly float _timeStep;
    private readonly BoundaryType _boundaryType;
    private float _mass;
    private bool _isActive;
    private Vector3 _position;
    private Vector3 _velocity;
    private float _radius;
    private readonly int _xDimension, _yDimension, _zDimension;

    public Vector3 Position { get => _position; }

    public ClassicalSystem(
        int x, int y, int z,
        BoundaryType boundaryType,
        float timeStep, float spaceStep,
        float mass)
    {
        _timeStep = timeStep;
        _xDimension = x;
        _yDimension = y;
        _zDimension = z;
        _boundaryType = boundaryType;
        _mass = mass;
        _isActive = true;
    }

    public void InitializeParticle(
        float xPosition, float yPosition, float zPosition,
        float xMomentumWavenumber, float yMomentumWavenumber, float zMomentumWavenumber,
        float radius)
    {
        _position = new Vector3(xPosition, yPosition, zPosition);
        _velocity = new Vector3(
            xMomentumWavenumber,
            yMomentumWavenumber,
            zMomentumWavenumber);
        _radius = radius;
    }

    public void HandleBoundaries()
    {
        if (!_isActive) return;

        switch (_boundaryType)
        {
            case BoundaryType.Reflective:
                ReflectiveBoundary();
                break;
            case BoundaryType.Absorbing:
                AbsorbingBoundary();
                break;
            case BoundaryType.Periodic:
                PeriodicBoundary();
                break;
        }
    }

    private void ReflectiveBoundary()
    {
        Vector3 newVelocity = _velocity;

        if (_position.X - _radius < 0 || _position.X + _radius > _xDimension) newVelocity.X = -newVelocity.X;
        if (_position.Y - _radius < 0 || _position.Y + _radius > _yDimension) newVelocity.Y = -newVelocity.Y;
        if (_position.Z - _radius < 0 || _position.Z + _radius > _zDimension) newVelocity.Z = -newVelocity.Z;

        _velocity = newVelocity;
    }

    private void AbsorbingBoundary()
    {
        if (_position.X < 0 || _position.X > _xDimension ||
            _position.Y < 0 || _position.Y > _yDimension ||
            _position.Z < 0 || _position.Z > _zDimension)
        {
            _isActive = false;
        }
    }

    private void PeriodicBoundary()
    {
        _position = new Vector3(
            (_position.X + _xDimension) % _xDimension,
            (_position.Y + _yDimension) % _yDimension,
            (_position.Z + _zDimension) % _zDimension);
    }

    public Vector3 CalculateForce(Vector3 position)
    {
        return Vector3.Zero;
    }

    public void Update()
    {
        if (!_isActive) return;

        Vector3 force = CalculateForce(_position);
        _velocity += force / _mass * _timeStep;
        _position += _velocity * _timeStep;

        HandleBoundaries();
    }
}