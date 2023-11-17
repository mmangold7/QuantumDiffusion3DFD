const cacheName = 'quantum-diffusion-cache-v1';
const assetsToCache = [
    './',  // Alias for 'index.html'
    'css/bootstrap/bootstrap.min.css',
    'css/app.css',
    'css/index.css',
    'favicon.png',
    'QuantumDiffusion3DFD.Client.styles.css',
    'manifest.json',
    'js/index.bundle.js',
    '_framework/blazor.webassembly.js',
    '_framework/wasm/mono.js',
    '_framework/wasm/dotnet.wasm',
    // Include additional resources required by Blazor
    // and any other assets like images, fonts, etc.
];

self.addEventListener('install', (event) => {
    // Pre-cache all required assets
    event.waitUntil(
        caches.open(cacheName)
            .then((cache) => {
                return cache.addAll(assetsToCache);
            })
    );
});

self.addEventListener('activate', (event) => {
    // Clean up old caches
    event.waitUntil(
        caches.keys().then((cacheNames) => {
            return Promise.all(
                cacheNames.filter((name) => name !== cacheName)
                    .map((name) => caches.delete(name))
            );
        })
    );
});

self.addEventListener('fetch', (event) => {
    // Serve cached content when offline
    event.respondWith(
        caches.match(event.request).then((response) => {
            return response || fetch(event.request);
        })
    );
});