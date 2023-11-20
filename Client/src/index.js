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

    let cubeIndex = 0;

    for (let x = 0; x < dimensions.x; x++) {
        for (let y = 0; y < dimensions.y; y++) {
            for (let z = 0; z < dimensions.z; z++) {
                const cubeSize = 10;
                const cubeGeometry = new THREE.BoxGeometry(cubeSize, cubeSize, cubeSize);
                const cubeMaterial = new THREE.MeshBasicMaterial({
                    color: "#000000",
                    transparent: true,
                    opacity: 0
                });

                const cube = new THREE.Mesh(cubeGeometry, cubeMaterial);
                cube.position.set(
                    x * spacing * spacingScaleFactor,
                    y * spacing * spacingScaleFactor,
                    z * spacing * spacingScaleFactor
                );

                cube.name = `cube-${cubeIndex}`; // Set a unique name for each cube
                scene.add(cube);
                cubeIndex++;
            }
        }
    }

    const wireFrameBox = new THREE.BoxGeometry(dimensions.x * singleDimensionScaleFactor, dimensions.y * singleDimensionScaleFactor, dimensions.z * singleDimensionScaleFactor);
    const boxEdges = new THREE.EdgesGeometry(wireFrameBox);
    const boxEdgeMaterial = new THREE.LineBasicMaterial({ color: 0xffffff, linewidth: 2 });
    const wireFrame = new THREE.LineSegments(boxEdges, boxEdgeMaterial);
    wireFrame.position.set(position.x, position.y, position.z);
    scene.add(wireFrame);

    animate();
}

function updateThreeJsScene(updatedData) {
    updatedData.forEach(data => {
        const { index, color, opacity } = data;

        let object = scene.getObjectByName(`cube-${index}`);
        if (object && object.isMesh) {
            object.material.color.set(color);
            object.material.opacity = opacity;
            object.material.transparent = true;
            object.material.blending = THREE.AdditiveBlending;
            object.material.needsUpdate = true;
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

async function loadWasmModule() {
    const wasmModule = await fetch('wasm/WASMkissFFT.wasm').then(response =>
        response.arrayBuffer()
    ).then(bytes =>
        WebAssembly.instantiate(bytes, {})
    );
    return wasmModule.instance.exports;
}

async function performFFT(data) {
    const wasmExports = await loadWasmModule();
    return wasmExports.fft(data);
}

export { initializeThreeJs, updatePointSize, updateThreeJsScene, performFFT };