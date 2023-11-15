var scene, camera, renderer, particleSystem;

function initializeThreeJs(dimensions, spacing) {
    const canvas = document.getElementById('threejs-canvas');
    if (!canvas) {
        console.error('Canvas element not found');
        return;
    }

    renderer = new THREE.WebGLRenderer({ canvas: canvas });

    const fov = 75;
    const aspect = 2;  // the canvas default
    const near = 0.1;
    const far = 1000;
    camera = new THREE.PerspectiveCamera(fov, aspect, near, far);
    camera.position.z = 150;

    scene = new THREE.Scene();

    const numParticles = dimensions.x * dimensions.y * dimensions.z;
    const positions = new Float32Array(numParticles * 3);
    const colors = new Float32Array(numParticles * 3);
    const color = new THREE.Color();

    let index = 0;
    for (let x = 0; x < dimensions.x; x++) {
        for (let y = 0; y < dimensions.y; y++) {
            for (let z = 0; z < dimensions.z; z++) {
                const i3 = index * 3;
                positions[i3 + 0] = x * spacing * 100; // x
                positions[i3 + 1] = y * spacing * 100; // y
                positions[i3 + 2] = z * spacing * 100; // z

                color.setRGB(1, 0, 0); // Default color (red)
                colors[i3 + 0] = color.r;
                colors[i3 + 1] = color.g;
                colors[i3 + 2] = color.b;

                index++;
            }
        }
    }

    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));

    const material = new THREE.PointsMaterial({
        size: 1.5,
        vertexColors: true
    });

    particleSystem = new THREE.Points(geometry, material);
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
    //console.log(probabilityData);

    const colors = particleSystem.geometry.attributes.color.array;

    // Normalize the probability data
    const maxProbability = Math.max(...probabilityData);

    for (let i = 0; i < probabilityData.length; i++) {
        const intensity = probabilityData[i] / maxProbability;
        const i3 = i * 3;
        colors[i3 + 0] = intensity; // Red channel
        colors[i3 + 1] = 0;         // Green channel
        colors[i3 + 2] = 0;         // Blue channel
    }

    particleSystem.geometry.attributes.position.needsUpdate = true;
    particleSystem.geometry.attributes.color.needsUpdate = true;
}
