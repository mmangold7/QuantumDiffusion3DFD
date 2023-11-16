var scene, camera, renderer, particleSystem;

function initializeThreeJs(dimensions, spacing) {
    const canvas = document.getElementById('threejs-canvas');
    if (!canvas) {
        console.error('Canvas element not found');
        return;
    }

    renderer = new THREE.WebGLRenderer({ canvas: canvas });

    const fov = 75;
    const aspect = canvas.clientWidth / canvas.clientHeight;
    const near = 0.1;
    const far = 1000;
    camera = new THREE.PerspectiveCamera(fov, aspect, near, far);
    camera.position.z = 150;

    scene = new THREE.Scene();

    const numParticles = dimensions.x * dimensions.y * dimensions.z;
    const positions = new Float32Array(numParticles * 3);
    const probabilities = new Float32Array(numParticles); // For the probabilities attribute

    let index = 0;
    for (let x = 0; x < dimensions.x; x++) {
        for (let y = 0; y < dimensions.y; y++) {
            for (let z = 0; z < dimensions.z; z++) {
                const i3 = index * 3;
                positions[i3 + 0] = x * spacing * 100; // x
                positions[i3 + 1] = y * spacing * 100; // y
                positions[i3 + 2] = z * spacing * 100; // z

                probabilities[index] = 1.0; // Initialize with a default value

                index++;
            }
        }
    }

    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute('probability', new THREE.BufferAttribute(probabilities, 1));

    const shaderMaterial = new THREE.ShaderMaterial({
        uniforms: {
            color: { value: new THREE.Color(0xff0000) },
            maxProbability: { value: 1.0 },
            probabilities: { value: [] },
            pointSize: { value: 10.0 } // Default point size
        },
        vertexShader: `
        attribute float probability;
        uniform float maxProbability;
        uniform float pointSize; // Use uniform for point size
        varying float vAlpha;
        void main() {
            vAlpha = probability / maxProbability;

            // Calculate distance from the camera
            vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
            float distance = length(mvPosition.xyz);

            // Adjust point size based on distance
            gl_PointSize = pointSize * vAlpha / distance * 300.0; // Adjust factor as needed

            gl_Position = projectionMatrix * mvPosition;
        }
    `,
        fragmentShader: `
        uniform vec3 color;
        varying float vAlpha;
        void main() {
            if (dot(gl_PointCoord - vec2(0.5), gl_PointCoord - vec2(0.5)) > 0.25) {
                discard;
            }
            gl_FragColor = vec4(color, vAlpha); // Entire point has uniform alpha based on probability
        }
    `,
        blending: THREE.AdditiveBlending,
        depthTest: false,
        transparent: true,
        vertexColors: true
    });

    particleSystem = new THREE.Points(geometry, shaderMaterial);
    scene.add(particleSystem);

    // Event listeners for mouse interaction
    let isLeftDragging = false;
    let isRightDragging = false;
    let isDragging = false;
    let previousMousePosition = {
        x: 0,
        y: 0
    };

    canvas.addEventListener('mousedown', (event) => {
        if (event.button === 0) {
            isLeftDragging = true;
        } else if (event.button === 2) {
            isRightDragging = true;
        }
        isDragging = isLeftDragging || isRightDragging;
        previousMousePosition.x = event.clientX;
        previousMousePosition.y = event.clientY;
    });

    canvas.addEventListener('mouseup', () => {
        isLeftDragging = false;
        isRightDragging = false;
        isDragging = false;
    });

    canvas.addEventListener('mousemove', (event) => {
        if (!isDragging) return;

        const deltaMousePosition = {
            x: event.clientX - previousMousePosition.x,
            y: event.clientY - previousMousePosition.y
        };

        const sensitivity = 0.002;

        if (isLeftDragging) {
            // Adjust rotation based on mouse movement
            camera.rotation.x += deltaMousePosition.y * sensitivity;
            camera.rotation.y += deltaMousePosition.x * sensitivity;
        } else if (isRightDragging) {
            // "Strafe" through the plane created by the current rotation and position
            const moveSpeed = 0.5;
            const moveDirection = new THREE.Vector3(-deltaMousePosition.x * moveSpeed, 0, -deltaMousePosition.y * moveSpeed);
            const rotationMatrix = new THREE.Matrix4().makeRotationFromEuler(camera.rotation);
            moveDirection.applyMatrix4(rotationMatrix);
            camera.position.add(moveDirection);
        }

        previousMousePosition.x = event.clientX;
        previousMousePosition.y = event.clientY;
    });

    // Scroll event listener
    canvas.addEventListener('wheel', (event) => {
        event.preventDefault(); // Prevent the default scrolling behavior
        const moveSpeed = 0.1; // Adjust this value for movement sensitivity
        const moveAmount = event.deltaY * moveSpeed;

        // Calculate the movement vector based on the camera's rotation
        const rotation = camera.rotation.clone();
        const movement = new THREE.Vector3(0, 0, moveAmount);
        movement.applyEuler(rotation);

        // Move the camera position
        camera.position.add(movement);
    });

    // Prevent the right-click context menu
    canvas.addEventListener('contextmenu', (event) => {
        event.preventDefault();
    });

    animate();
}

function updatePointSize(newPointSize) {
    if (particleSystem && particleSystem.material && particleSystem.material.uniforms.pointSize) {
        particleSystem.material.uniforms.pointSize.value = newPointSize;
        particleSystem.material.needsUpdate = true;
    }
}

function animate() {
    requestAnimationFrame(animate);
    render();
}

function render() {
    if (resizeRendererToDisplaySize(renderer)) {
        const canvas = renderer.domElement;
        camera.aspect = canvas.clientWidth / canvas.clientHeight;
        camera.updateProjectionMatrix();
    }
    renderer.render(scene, camera);
}

function resizeRendererToDisplaySize(renderer) {
    const canvas = renderer.domElement;
    const width = canvas.clientWidth;
    const height = canvas.clientHeight;
    const needResize = canvas.width !== width || canvas.height !== height;
    if (needResize) {
        renderer.setSize(width, height, false);
    }
    return needResize;
}

function updateThreeJsScene(probabilityData) {
    if (particleSystem) {
        const probabilities = particleSystem.geometry.attributes.probability.array;
        let maxProbability = Math.max(...probabilityData);

        // If your maxProbability is zero (which can't be used for division), then set it to 1
        if (maxProbability === 0) maxProbability = 1;

        for (let i = 0; i < probabilityData.length; i++) {
            probabilities[i] = probabilityData[i] / maxProbability; // Normalized probability
        }

        // Update the maximum probability uniform
        particleSystem.material.uniforms.maxProbability.value = maxProbability;

        // Mark the probabilities attribute as needing an update
        particleSystem.geometry.attributes.probability.needsUpdate = true;
    }
}
