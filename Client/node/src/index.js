import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls';

var scene, camera, renderer, particleSystem;
var controls;

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
    camera.position.z = 150;

    controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.1;
    controls.rotateSpeed = 0.1;
    controls.zoomSpeed = 1.0;
    controls.panSpeed = 0.8;

    scene = new THREE.Scene();

    const numParticles = dimensions.x * dimensions.y * dimensions.z;
    const positions = new Float32Array(numParticles * 3);
    const probabilities = new Float32Array(numParticles);

    let index = 0;
    for (let x = 0; x < dimensions.x; x++) {
        for (let y = 0; y < dimensions.y; y++) {
            for (let z = 0; z < dimensions.z; z++) {
                const i3 = index * 3;
                positions[i3 + 0] = x * spacing * spacingScaleFactor; // x
                positions[i3 + 1] = y * spacing * spacingScaleFactor; // y
                positions[i3 + 2] = z * spacing * spacingScaleFactor; // z

                probabilities[index] = 1.0;

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
            pointSize: { value: 10.0 }
        },
        vertexShader: `
        attribute float probability;
        uniform float maxProbability;
        uniform float pointSize; // Use uniform for point size
        varying float vAlpha;
        varying float vProbability; // Declare vProbability as varying

        void main() {
            vAlpha = probability / maxProbability;
            vProbability = probability; // Pass the probability to the fragment shader

            // Calculate distance from the camera
            vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
            float distance = length(mvPosition.xyz);

            // Adjust point size based on distance
            gl_PointSize = pointSize * vAlpha / distance * 300.0; // Adjust factor as needed

            gl_Position = projectionMatrix * mvPosition;
        }
    `,
        fragmentShader: `
        varying float vProbability;

        vec3 interpolateColor(float value) {
            vec3 colorStops[7] = vec3[](vec3(1, 0, 0), vec3(1, 0.5, 0), vec3(1, 1, 0), vec3(0, 1, 0), vec3(0, 0, 1), vec3(0.29, 0, 0.51), vec3(0.56, 0, 1)); // ROYGBIV
            float scaledValue = value * 6.0;
            int index = int(scaledValue);
            vec3 color1 = colorStops[index % 7];
            vec3 color2 = colorStops[(index + 1) % 7];
            float frac = fract(scaledValue);
            return mix(color1, color2, frac);
        }

        void main() {
            if (dot(gl_PointCoord - vec2(0.5), gl_PointCoord - vec2(0.5)) > 0.25) {
                discard;
            }
            
            vec3 color = interpolateColor(vProbability);
            gl_FragColor = vec4(color, 1.0);
        }
    `,
        blending: THREE.AdditiveBlending,
        depthTest: false,
        transparent: true,
        vertexColors: true
    });

    particleSystem = new THREE.Points(geometry, shaderMaterial);
    scene.add(particleSystem);

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

function updateThreeJsScene(probabilityData) {
    if (particleSystem) {
        const probabilities = particleSystem.geometry.attributes.probability.array;
        let maxProbability = Math.max(...probabilityData);

        // Zero can't be used for division
        if (maxProbability === 0) maxProbability = 1;

        for (let i = 0; i < probabilityData.length; i++) {
            probabilities[i] = probabilityData[i] / maxProbability; // Normalized probability
        }

        particleSystem.material.uniforms.maxProbability.value = maxProbability;
        particleSystem.geometry.attributes.probability.needsUpdate = true;
    }
}

function updatePointSize(newPointSize) {
    if (particleSystem && particleSystem.material && particleSystem.material.uniforms.pointSize) {
        particleSystem.material.uniforms.pointSize.value = newPointSize;
        particleSystem.material.needsUpdate = true;
    }
}

//hacky workaround until I figure out why webpack isn't putting these on the window like it's configured to
window.QuantumInterop = {
    initializeThreeJs,
    updatePointSize,
    updateThreeJsScene
};

export { initializeThreeJs, updatePointSize, updateThreeJsScene };