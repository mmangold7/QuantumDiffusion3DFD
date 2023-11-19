import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls';

var scene, camera, renderer, particleSystem, controls;

function initializeThreeJs(dimensions, spacing) {
    const canvas = document.getElementById('threejs-canvas');
    var spacingScaleFactor = 100;
    var singleDimensionScaleFactor = spacingScaleFactor / 10;
    var position = {
        x: (dimensions.x - 1) * singleDimensionScaleFactor / 2,
        y: (dimensions.y - 1) * singleDimensionScaleFactor / 2,
        z: (dimensions.z - 1) * singleDimensionScaleFactor / 2
    }
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
    controls = new OrbitControls(camera, renderer.domElement);
    camera.position.x = position.x;
    camera.position.x = position.y;
    camera.position.x = position.z + 100;
    controls.update();
    controls.target.set(position.x, position.y, position.z);

    controls.enableDamping = true;
    controls.dampingFactor = 0.1;
    controls.rotateSpeed = 0.1;
    controls.zoomSpeed = 1.0;
    controls.panSpeed = 0.8;

    scene = new THREE.Scene();

    for (let x = 0; x < dimensions.x; x++) {
        for (let y = 0; y < dimensions.y; y++) {
            for (let z = 0; z < dimensions.z; z++) {
                const probability = 0.5;

                const cubeSize = 10;
                const cubeGeometry = new THREE.BoxGeometry(cubeSize, cubeSize, cubeSize);
                const color = interpolateColor(probability);
                const cubeMaterial = new THREE.MeshBasicMaterial({
                    color: color,
                    transparent: true,
                    opacity: probability // Opacity based on probability
                });

                const cube = new THREE.Mesh(cubeGeometry, cubeMaterial);
                cube.position.set(
                    x * spacing * spacingScaleFactor,
                    y * spacing * spacingScaleFactor,
                    z * spacing * spacingScaleFactor
                );

                scene.add(cube);
            }
        }
    }

    const boxGeometry = new THREE.BoxGeometry(dimensions.x * singleDimensionScaleFactor, dimensions.y * singleDimensionScaleFactor, dimensions.z * singleDimensionScaleFactor);
    const edgesGeometry = new THREE.EdgesGeometry(boxGeometry);
    const lineMaterial = new THREE.LineBasicMaterial({ color: 0xffffff, linewidth: 2 });
    const wireframe = new THREE.LineSegments(edgesGeometry, lineMaterial);
    wireframe.position.set(position.x, position.y, position.z);
    scene.add(wireframe);

    //const arrowGeometry = new THREE.ConeGeometry(0.5, 1, 32); // Small cone to represent an arrow
    //const arrowMaterial = new THREE.MeshBasicMaterial({ color: 0x00ff00 }); // Green color
    //const arrow = new THREE.Mesh(arrowGeometry, arrowMaterial);
    //arrow.rotation.x = Math.PI / 2; // Rotate to point upwards
    //arrow.visible = false; // Initially hidden
    //scene.add(arrow);

    animate();
}

function updateThreeJsScene(probabilityData) {
    let maxProbability = Math.max(...probabilityData);
    if (maxProbability === 0) maxProbability = 1; // Avoid division by zero

    let cubeIndex = 0;

    scene.traverse(function (object) {
        if (object.isMesh) {
            const probability = probabilityData[cubeIndex] / maxProbability;
            const color = interpolateColor(probability);
            const opacity = probabilityToOpacity(probability);

            object.material.color.set(color);
            object.material.opacity = opacity;
            object.material.transparent = true;
            object.material.blending = THREE.AdditiveBlending;
            object.material.needsUpdate = true;

            cubeIndex++;
        }
    });
}

function updatePointSize(newPointSize) {
    if (particleSystem && particleSystem.material && particleSystem.material.uniforms.pointSize) {
        particleSystem.material.uniforms.pointSize.value = newPointSize;
        particleSystem.material.needsUpdate = true;
    }
}

function animate() {
    requestAnimationFrame(animate);
    controls.update();
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

function probabilityToOpacity(probability) {
    // Adjust this function to fine-tune how opacity scales with probability
    //let opacity = probability; // Direct linear mapping
    //let opacity = Math.sqrt(probability); // Square root mapping
    //let opacity = probability * probability; // Squaring the probability
    //let opacity = Math.log(probability + 1) / Math.log(2); // Log base 2
    let opacity = 1 / (1 + Math.exp(-10 * (probability - 0.5))); // Sigmoid function
    return opacity;
}

function interpolateColor(value) {
    // Function to return a color based on the value
    // Modify this function to get the desired color mapping
    const hue = (1 - value) * 240; // From blue (low) to red (high)
    return `hsl(${hue}, 100%, 50%)`;
}

function interpolateColorRainbow(value) {
    const colorStops = [
        new THREE.Color(1, 0, 0), // Red
        new THREE.Color(1, 0.5, 0), // Orange
        new THREE.Color(1, 1, 0), // Yellow
        new THREE.Color(0, 1, 0), // Green
        new THREE.Color(0, 0, 1), // Blue
        new THREE.Color(0.29, 0, 0.51), // Indigo
        new THREE.Color(0.56, 0, 1) // Violet
    ];

    let scaledValue = value * (colorStops.length - 1);
    let index = Math.floor(scaledValue);
    let frac = scaledValue - index;

    let color1 = colorStops[index];
    let color2 = colorStops[Math.min(index + 1, colorStops.length - 1)];

    return color1.clone().lerp(color2, frac);
}

export { initializeThreeJs, updatePointSize, updateThreeJsScene };