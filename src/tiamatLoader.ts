import * as THREE from 'three'
import * as base from "./libs/base";
import * as utils from "./libs/utils";

const DNA_BASES = new Map<string,number>([
    ['A', 0],
    ['C', 1],
    ['G', 2],
    ['T', 3],
    ['U', 3]]
);

function loadTiamat(source_file: string, tiamat_version: number, isDNA=true, default_val="R") {
    let tiamat_version_fudge: number;
    if (tiamat_version === 1) {
        if (isDNA) {
            tiamat_version_fudge = 1.2;
        } else { 
            tiamat_version_fudge = 1.6;
        }
    } else if (tiamat_version === 2) {
        tiamat_version_fudge = 1
    } else {
        throw `tiamat_version should be either '1' or '2' (got '${tiamat_version}' instead`;
    }

    const data = JSON.parse(source_file);

    // Extract and initialize the bases
    const bases = data['bases'].map(b=>new TiamatBase(b, isDNA, tiamat_version_fudge, default_val));

    // get a list of the strands in the system
    const system = create_system(bases, isDNA);

    return system;
}


// creates a list of all strands in the respective tiamat file 
function create_system(bases: TiamatBase[], isDNA=true) {
    // lookup bases by id
    let bases_by_id = new Map(bases.map(b=>[b.id, b]));
    
    // Calc base vectors
    bases.forEach(b => {
        b.set_base_config(bases_by_id);
    });

    let box_size = 1;
    // Make sure box fits bases
    for (const b of bases) {
        for(const v of b.cm_pos.toArray())  {
            if (Math.abs(v) > box_size/2) {
                box_size = Math.abs(v)*2
            }
        }
    }
    box_size *= 1.5; // Add some margin
    box_size = Math.round(box_size);
    
    // Initialize system
    let sys = new base.System(new THREE.Vector3(box_size, box_size, box_size));
    sys.isDNA = isDNA;

    let included = new Set<number>();

    //  5' ends don't have bases downstream
    let fivePrimes = bases.filter(b => b.obj['down'] === null);

    for (const start_base of fivePrimes) {
        let strand = new base.Strand();
        let b = start_base;
        // follows every strand up the backbone (5'->3')
        while (b !== undefined) {
            let nuc = new base.Nucleotide(b.cm_pos, b.a1, b.a3, b.val);
            b.nuc = nuc;
            strand.add_nucleotide(nuc);
            included.add(b.id);
            b = b.getNeighbor(bases_by_id, 'up');
        }
        sys.add_strand(strand);
    }

    let extra_strands = get_circular_strands(bases, bases_by_id, included);
    sys.add_strands(extra_strands);

    // Set basepairs
    bases.forEach(b=>{
        if (b.obj['across'] !== null) {
            b.nuc.pair = bases_by_id.get(b.obj['across']).nuc;
            // If the tiamat file did not specify any type,
            // set it to be complementary to the paired type
            if (b.obj['type'] === undefined) {
                b.nuc._base = 3 - b.nuc.pair._base;
            }
        }
    });
    base.Logger.log(`Added base pairs strands`);

    sys.calc_clusters();

    return sys;
}

function get_circular_strands(bases: TiamatBase[], bases_by_id: Map<number, TiamatBase>, already_included: Set<number>) {
    let cs = bases.filter(b=>!already_included.has(b.id));
    let strands = []

    while (cs.length > 0) {
        let strand = new base.Strand();
        let base_lst_ids = new Set();
        let b = cs[0];
        // Loop until we come full circle
        // (and try to add the same base again)
        while (!base_lst_ids.has(b.id)) {
            let nuc = new base.Nucleotide(b.cm_pos, b.a1, b.a3, b.val);
            b.nuc = nuc;
            strand.add_nucleotide(nuc);
            base_lst_ids.add(b.id);
            b = b.getNeighbor(bases_by_id, 'up');
        }
        strands.push(strand);
        cs = cs.filter(b=>!base_lst_ids.has(b.id));
    }
    return strands
}

class TiamatBase {
    static scale = 1 / 0.85;

    id: number;
    up: TiamatBase;
    down: TiamatBase;
    across: TiamatBase;
    val: number;
    tiamat_pos: THREE.Vector3;
    cm_pos: THREE.Vector3;
    a1: THREE.Vector3;
    a3: THREE.Vector3;
    obj: Object;
    tiamat_version_fudge: number;
    isDNA: boolean;
    nuc: base.Nucleotide;

    constructor(obj: Object, isDNA=true, tiamat_version_fudge: number, default_val="R") {
        this.id = obj['id'];
        this.isDNA = isDNA;
        this.tiamat_version_fudge = tiamat_version_fudge;
        try {
            this.val = DNA_BASES.get(obj['type'][0]);
        } catch {
            if (default_val == 'R') {
                this.val = utils.randint(0,4) as number;
            } else {
                this.val = DNA_BASES.get(default_val);
                console.assert(this.val !== undefined, `Default base ${default_val} not supported, use one of ${DNA_BASES}`);
            }
            const baseNames = ['A', 'C', 'G', isDNA ? 'T':'U']
            base.Logger.log(`WARNING: base n.${this.id} has no associated type, setting it to '${baseNames[this.val]}'`);
        }
        this.tiamat_pos = new THREE.Vector3().fromArray(obj['position']);
        this.tiamat_pos.multiplyScalar(TiamatBase.scale);
        this.obj = obj;
    }

    getNeighborId(relation: string): number {
        let id = this.obj[relation];
        if (id === undefined) {
            id = -1;
        }
        return id;
    }

    getNeighbor(bases_by_id: Map<number, TiamatBase>, ...relations: string[]) {
        let b = this as TiamatBase;
        for (const relation of relations) {
            b = bases_by_id.get(b.obj[relation]);
            if (b === undefined) {
                break;
            }
        }
        return b;
    }

    hasNeighbor(bases_by_id: Map<number, TiamatBase>, ...relations: string[]) {
        return this.getNeighbor(bases_by_id, ...relations) !== undefined;
    }

    getNeighborPos(bases_by_id: Map<number, TiamatBase>, ...relations: string[]) {
        let b = this.getNeighbor(bases_by_id, ...relations);
        let pos: THREE.Vector3;
        if (b === undefined) {
            pos = new THREE.Vector3(0, 0, 1).multiplyScalar(TiamatBase.scale);
        } else {
            pos = b.tiamat_pos;
        }
        return pos;
    }

    // Call after all bases have been initialized
    set_base_config(bases_by_id: Map<number, TiamatBase>) {
        let base_vector = this.tiamat_pos;

        let paring_base_vector = this.getNeighborPos(bases_by_id, 'across');
        let up_base_vector = this.getNeighborPos(bases_by_id, 'up');
        let down_base_vector = this.getNeighborPos(bases_by_id, 'down');
        let paring_base3_vector = this.getNeighborPos(bases_by_id, 'up', 'across');
        let paring_base5_vector = this.getNeighborPos(bases_by_id, 'down', 'across');

        // three backbone vectors, A-base_vector, B-paring_base_vector,
        // A5up_base_vector, A3down_base_vector
        let backbone_A_to_backbone_B = paring_base_vector.clone().sub(base_vector).normalize();
        let backbone_A_to_posA5_neighbor = up_base_vector.clone().sub(base_vector).normalize();
        let backbone_A_to_posA3_neighbor = down_base_vector.clone().sub(base_vector).normalize();
        let backboneA_to_backbone_B3 = paring_base3_vector.clone().sub(base_vector).normalize();
        let backboneA_to_backbone_B5 = paring_base5_vector.clone().sub(base_vector).normalize();

        // do we have a duplex across the 3' end
        if (this.hasNeighbor(bases_by_id, 'down', 'across')) {
            if (this.isDNA) {
                [this.a1, this.a3, this.cm_pos] = neighbor3_cal_vector(backbone_A_to_backbone_B, backbone_A_to_posA3_neighbor, backboneA_to_backbone_B5)
            } else {
                [this.a1, this.a3, this.cm_pos] = neighbor3_cal_vector_RNA(backbone_A_to_backbone_B, backbone_A_to_posA3_neighbor, backboneA_to_backbone_B5)
            }
            this.cm_pos.add(base_vector.clone().multiplyScalar(this.tiamat_version_fudge));
        // do we have a duplex across the 5' end
        } else if (this.hasNeighbor(bases_by_id, 'up', 'across')) {
            if (this.isDNA) {
                [this.a1, this.a3, this.cm_pos] = neighbor5_cal_vector(backbone_A_to_backbone_B, backbone_A_to_posA5_neighbor, backboneA_to_backbone_B3)
            } else {
                [this.a1, this.a3, this.cm_pos] = neighbor5_cal_vector_RNA(backbone_A_to_backbone_B, backbone_A_to_posA5_neighbor, backboneA_to_backbone_B3);
            }
            this.cm_pos.add(base_vector.clone().multiplyScalar(this.tiamat_version_fudge));
        } else {
            // single stranded case, not really treated (could randomize orientations?)
            this.cm_pos = base_vector
            this.a3 = new THREE.Vector3(0., 0., 1.);
            this.a1 = new THREE.Vector3(0., 1., 0.);
        }
    }
}

function calcVec(a: THREE.Vector3, b: THREE.Vector3, c: THREE.Vector3,
    af1: number, bf1: number, cf1: number,
    af2: number, bf2: number, cf2: number
    ) {
    let calc = (cA: number, cB: number, cC: number) => 
        a.clone().multiplyScalar(cA).add(
        b.clone().multiplyScalar(cB).add(
        c.clone().multiplyScalar(cC))
    ).normalize();
    const cm = calc(af1, bf1, cf1);
    const a3 = calc(af2, bf2, cf2);
    const a1 = a.clone();
    return [a1, a3, cm];
}

// used to calculate the center of mass 
function neighbor3_cal_vector(AB: THREE.Vector3, AA3: THREE.Vector3, AB5: THREE.Vector3) {
    return calcVec(AB, AA3, AB5,
        -0.13079674, -0.22543211, 0.62949112,
        2.69498211, -1.04531113, -2.30531223
    );
}

// used to calculate the center of mass 
function neighbor5_cal_vector(AB: THREE.Vector3, AA5: THREE.Vector3, AB3: THREE.Vector3) {
    return calcVec(AB, AA5, AB3, 
        0.81079674, 0.22543211, - 0.50262804,
        -2.12846367, 0.82557385, 2.33064701
    );
}

function neighbor3_cal_vector_RNA(AB: THREE.Vector3, AA3: THREE.Vector3, AB5: THREE.Vector3) {
    return calcVec(AB, AA3, AB5,
        -0.28102082, -0.25891019, 0.84990909,
        2.34763359, -1.1627428, -1.63537381
    );
}
function neighbor5_cal_vector_RNA(AB: THREE.Vector3, AA5: THREE.Vector3, AB3: THREE.Vector3) {
    return calcVec(AB, AA5, AB3,
        0.85540635, 0.30569283, -0.44567833,
        -1.60523423, 0.58820649, 2.00150202
    );
}

export {loadTiamat}