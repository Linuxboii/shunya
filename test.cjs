// node test.cjs — checks the zero finder against published values and counts
const assert = require('node:assert');
const Z = require('./zeta.js');
const known = [14.134725142, 21.022039639, 25.010857580, 30.424876126, 32.935061588, 37.586178159, 40.918719012, 43.327073281, 48.005150881, 49.773832478];
const found = Z.zeros(10, 50);
assert.equal(found.length, known.length);
found.forEach((z, i) => assert.ok(Math.abs(z - known[i]) < 1e-8, `zero ${i + 1}: ${z} vs ${known[i]}`));
assert.equal(Z.zeros(0, 1000).length, 649, 'N(1000)');
assert.equal(Z.zeros(0, 10000).length, 10142, 'N(10000)');
const lehmer = Z.zeros(7005, 7005.2);
assert.equal(lehmer.length, 2, "Lehmer's close pair near t = 7005.1");
console.log('ok: first 10 zeros to 1e-8, N(1000) = 649, N(10000) = 10142, Lehmer pair resolved');
