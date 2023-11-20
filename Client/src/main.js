﻿import { initializeThreeJs, updatePointSize, updateThreeJsScene, performFFT } from './index';

window.addEventListener('load', () => {
    if ('serviceWorker' in navigator) {
        const swFile = process.env.ENVIRONMENT === 'production' ? 'service-worker.prod.js' : 'service-worker.dev.js';
        navigator.serviceWorker.register(swFile);
    }
});

window.QuantumInterop = {
    initializeThreeJs,
    updatePointSize,
    updateThreeJsScene,
    performFFT
};