import {FLT_EPSILON} from "./base"
import * as THREE from 'three'

function sum(...values: number[]) {
    return values.reduce((a, b) => a + b, 0)
}

function range(low: number, high?: number, step=1) {
    if (high === undefined) {
        high = low;
        low = 0;
    }
    let r = [];
    for (let i=low; i<high; i+=step) {
        r.push(i);
    }
    return r;
}

// Adepted from https://github.com/nvie/itertools.js
function* permutations<T>(iterable: Iterable<T>, r?: number): Iterable<Array<T>> {
    let pool = Array.from(iterable);
    let n = pool.length;
    let x = r === undefined ? n : r;

    if (x > n) {
        return;
    }

    let indices: Array<number> = Array.from(range(n));
    let cycles: Array<number> = Array.from(range(n, n - x, -1));
    let poolgetter = (i) => pool[i];

    yield indices.slice(0, x).map(poolgetter);

    while (n > 0) {
        let cleanExit: boolean = true;
        for (let i of range(x - 1, -1, -1)) {
            cycles[i] -= 1;
            if (cycles[i] === 0) {
                indices = indices
                    .slice(0, i)
                    .concat(indices.slice(i + 1))
                    .concat(indices.slice(i, i + 1));
                cycles[i] = n - i;
            } else {
                let j: number = cycles[i];

                let [p, q] = [indices[indices.length - j], indices[i]];
                indices[i] = p;
                indices[indices.length - j] = q;
                yield indices.slice(0, x).map(poolgetter);
                cleanExit = false;
                break;
            }
        }

        if (cleanExit) {
            return;
        }
    }
}

// Behave like numpy.random.randint
function randint(low: number, high?: number, size=1): number | number[] {
    if (high === undefined) {
        high = low;
        low = 0;
    }
    var a = [];
    for (let i=0;i<size;i++){
        a[i] = low + Math.floor(Math.random()*(high-low));
    };
    if (size === 1) {
        return a[0];
    }
    return a;
}

function baseFromIUPAC(s: string, rna=false) {
    const T = rna ? 'T' : 'U';
    const r = (l: any[]) => l[Math.floor(Math.random()*l.length)];
    return [...s].map(c=>{
        if (['A', 'T', 'U', 'C', 'G'].includes(c)) {
            return c
        }
        switch (c.toUpperCase()) {
            case 'R': return r(['A', 'G']);         // Purine
            case 'Y': return r(['C',  T ]);         // Pyrimindine
            case 'M': return r(['A', 'C']);
            case 'K': return r(['G',  T ]);
            case 'S': return r(['C', 'G']);         // Strong
            case 'W': return r(['A',  T ]);         // Week
            case 'H': return r(['A', 'C',  T ]);    // Not G
            case 'B': return r(['C', 'G',  T ]);    // Not A
            case 'V': return r(['A', 'C', 'G']);    // Not U/T
            case 'D': return r(['A', 'G',  T ]);    // Not C
            case 'N': return r(['A', 'C', 'G', T]); // Ambiguous
            default:
                throw new Error(`Unknown base code ${c}`);
                break;
        }
    }).join('')
}

//From https://stackoverflow.com/a/16436975
function arraysEqual(a: any[], b: any[]) {
    if (a === b) return true;
    if (a === null || b === null) return false;
    if (a.length !== b.length) return false;
  
    for (var i = 0; i < a.length; ++i) {
      if (a[i] !== b[i]) return false;
    }
    return true;
  }

function get_angle(a, b) {
    /*
    Get angle between a,b

    >>> a =[0, 1, 0];
    >>> b =[0, 0, 1];
    >>> round(get_angle(a,b),3)
    1.571

    */
    let ab = a.dot(b);
    if (ab > (1. - FLT_EPSILON)) { 
        return 0;
    } else if (ab < (-1. + FLT_EPSILON)) { 
        return Math.PI;
    } else { 
        return Math.acos(ab)
    }
}


function get_orthonormalized_base(v1: THREE.Vector3, v2: THREE.Vector3, v3: THREE.Vector3): [THREE.Vector3, THREE.Vector3, THREE.Vector3] {
    v1 = v1.clone();
    v2 = v2.clone();
    v3 = v3.clone();

    const v1_norm2 = v1.dot(v1);
    const v2_v1 = v2.dot(v1);

    v2.sub(v1.clone().multiplyScalar((v2_v1 / v1_norm2)));

    const v3_v1 = v3.dot(v1);
    const v3_v2 = v3.dot(v2);
    const v2_norm2 = v2.dot(v2);

    v3.sub(v1.clone().multiplyScalar((v3_v1 / v1_norm2)).add(v2.clone().multiplyScalar((v3_v2 / v2_norm2))));

    v1.divideScalar(v1_norm2);
    v2.divideScalar(v2_norm2);
    v3.divideScalar(Math.sqrt(v3.dot(v3)));

    return [v1, v2, v3];
}


function get_random_vector_in_sphere(r=1) {
    const r2 = r*r;

    let rnd = ()=>(2*r*Math.random()-r);

    let v = new THREE.Vector3(rnd(), rnd(), rnd());
    while (v.dot(v) > r2) {
        v = new THREE.Vector3(rnd(), rnd(), rnd());
    }

    return v;
}


function get_random_vector() {
    let ransq = 1.;

    let ran1: number, ran2: number;
    while (ransq >= 1.) {
        ran1 = 1. - 2. * Math.random();
        ran2 = 1. - 2. * Math.random();
        ransq = ran1 * ran1 + ran2 * ran2;
    }

    const ranh = 2. * Math.sqrt(1. - ransq);
    return new THREE.Vector3(ran1 * ranh, ran2 * ranh, 1. - 2. * ransq);
}


function get_random_rotation_matrix() {
    let [v1, v2, v3] = get_orthonormalized_base(get_random_vector(), get_random_vector(), get_random_vector());

    let R = new THREE.Matrix3().set(
        v1.x, v1.y, v1.z,
        v2.x, v2.y, v2.z,
        v3.x, v3.y, v3.z,
    );
    //  rotations have det === 1
    if (R.determinant() < 0) {
        R = new THREE.Matrix3().set(
            v2.x, v2.y, v2.z,
            v1.x, v1.y, v1.z,
            v3.x, v3.y, v3.z
        );
    }

    return R;
}

export {sum, range, permutations, randint, baseFromIUPAC, arraysEqual, get_angle, get_orthonormalized_base, get_random_vector_in_sphere, get_random_vector, get_random_rotation_matrix};
