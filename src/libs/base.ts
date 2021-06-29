import * as THREE from 'three'
import * as utils from './utils'
/*
Utility functions.
base.py includes the classes: System, Strand, Nucleotide
    - Make initial configurations (generate.py)
    - Get detailed energy information (process_data/)
    - If you want to use it with oxRNA, you have to set environment variable OXRNA to 1  (export OXRNA=1) 
*/

const number_to_base = {0 : 'A', 1 : 'G', 2 : 'C', 3 : 'T'};

const base_to_number = {
    'A' : 0, 'a' : 0, 'G' : 1, 'g' : 1,
    'C' : 2, 'c' : 2, 'T' : 3, 't' : 3,
    'U' : 3, 'u' : 3, 'D' : 4
};

let FLT_EPSILON: number;
try {
    FLT_EPSILON = Number.EPSILON;
} catch {
    FLT_EPSILON = 2.2204460492503131e-16;
}
    
// oxDNA constants
const POS_BACK = -0.4;
const POS_STACK = 0.34;
const POS_BASE = 0.4;
const CM_CENTER_DS = POS_BASE + 0.2;
const FENE_R0_OXDNA = 0.7525;
const FENE_EPS = 2.0;

const LENGTH_FACT = 8.518;
const BASE_BASE = 0.3897628551303122;

const CREPY_COLOR_TABLE = ['red', 'blue', '0,0.502,0', '1,0.8,0', '0.2,0.8,1'];

// INT_POT_TOTAL = 0;
const INT_HYDR = 4;
const INT_STACK = 2;
const INT_CROSS_STACK = 5;
const INT_COAX_STACK = 6;
const INT_FENE = 0;
const INT_EXC_BONDED = 1;
const INT_EXC_NONBONDED = 3;

const H_CUTOFF = -0.1;

const MM_GROOVING = false;


class Logger {
    static DEBUG = 0;
    static INFO = 1;
    static WARNING = 2;
    static CRITICAL = 3;
    static debug_level = Logger.INFO;
    static messages = ["DEBUG", "INFO", "WARNING", "CRITICAL"];

    static logFunction: (msg: string) => void = console.log;

    static log(msg: string, level=Logger.INFO, additional?: string) {
        if (level === undefined) {
            level = Logger.INFO
          }
        if (level < Logger.debug_level) {
            return;
        }
        if (additional != undefined && Logger.debug_level === Logger.DEBUG) {
            Logger.logFunction(`${Logger.messages[level]}: ${msg} (additional info: '${additional}')`);
        }
        else {
            Logger.logFunction(`${Logger.messages[level]}: ${msg}`);
        }
    }

    static die(msg) {
        Logger.log(msg, Logger.CRITICAL)
        throw msg;
    }
}


class Nucleotide {
    /*
    Nucleotides compose Strands

    cm_pos --- Center of mass position
        Ex: [0, 0, 0]

    a1 --- Unit vector indicating orientation of backbone with respect to base
        Ex: [1, 0, 0]

    a3 --- Unit vector indicating orientation (tilting )of base with respect to backbone
        Ex: [0, 0, 1]

    base --- Identity of base, which must be designated with either numbers or
        letters (this is called type in the c++ code). Confusingly enough, this
        is similar to Particle.btype in oxDNA.
        
        Number: {0,1,2,3} or any number in between (-inf,-7) and (10, inf)
        To use specific sequences (or an alphabet large than four) one should
        start from the complementary pair 10 and -7. Complementary pairs are
        such that base_1 + base_2 = 3;;
        
        Letter: {A,G,T,C} (these will be translated to {0, 1, 3, 2}).
        
        These are set in the dictionaries: number_to_base, base_to_number

    pair --- Base-paired Nucleotide, used in oxView output
    cluster --- Cluster ID, Number, used in oxView output
    color --- Custom color, Number representing hex value, used in oxView output

    */

    static index = 0;
    index: number;
    cm_pos: THREE.Vector3;
    cm_pos_box: THREE.Vector3;
    _a1: THREE.Vector3;
    _a3: THREE.Vector3;
    _base: number;
    _L;
    _v;
    n3: number;
    next: number;
    pair: Nucleotide;
    cluster: number
    color: number;
    strand: number;

    constructor(cm_pos, a1, a3, base: number | string, 
        v=new THREE.Vector3(0., 0., 0.), L=new THREE.Vector3(0., 0., 0.),
        n3=-1,pair=undefined, cluster=undefined,color=undefined
    ) {
        this.index = Nucleotide.index;
        Nucleotide.index++
        this.cm_pos = cm_pos.clone();
        this._a1 = a1.clone().normalize();
        this._a3 = a3.clone().normalize();

        // Set random base if undefined
        if (base === undefined) {
            base = utils.randint(0,4) as number;
        }
        //  base should be an integer
        if (typeof base === "string") {
            try {
                base = base_to_number[base];
            } catch (e) {
                Logger.log(`Invalid base (${e})`);
            }
        }
        this._base = base as number;
        this._L = L;
        this._v = v;
        this.n3 = n3;
        this.next = -1;
        this.pair = pair;
        this.cluster = cluster;
        this.color = color;
    }

    /**
     * Get position of the base centroid
     * Note that cm_pos is the centrod of the backbone and base.
     * @returns the position of the base centroid
     */
    get pos_base(): THREE.Vector3 {
        return this.cm_pos.clone().add(this._a1.clone().multiplyScalar(POS_BASE));
    }

    get pos_stack (): THREE.Vector3 {
        return this.cm_pos.clone().add(this._a1.clone().multiplyScalar(POS_STACK));
    }

    get pos_back (): THREE.Vector3 {
        return this.cm_pos.clone().add(this._a1.clone().multiplyScalar(POS_BACK));
    }

    /**
     * Get the position of the backbone centroid relative to the centre of mass
        i.e. it will be a vector pointing from the c.o.m. to the backbone
     * @returns position of the backbone centroid relative to the centre of mass
     */
    get pos_back_rel(): THREE.Vector3 {
        return this.pos_back.clone().sub(this.cm_pos);
    }

    get a2 () {
        return this._a3.clone().cross(this._a1);
    }

    copy(disp?: THREE.Vector3, rot?: THREE.Matrix3) {
        let copy = new Nucleotide(
            this.cm_pos.clone(), this._a1.clone(), this._a3.clone(),
            this._base, this._L, this._v, this.n3, this.pair,
            this.cluster, this.color
        );
        if (disp !== undefined) {
            copy.translate(disp);
        }
        if (rot !== undefined) {
            copy.rotate(rot);
        }
        return copy;
    }

    translate(disp: THREE.Vector3) {
        this.cm_pos.add(disp);
        this.cm_pos_box.add(disp);
    }

    rotate(R: THREE.Matrix3, origin?: THREE.Vector3) {
        if (origin === undefined) {
            origin = this.cm_pos.clone();
        } 

        this.cm_pos.sub(origin).applyMatrix3(R).add(origin);
        this._a1.applyMatrix3(R);
        this._a3.applyMatrix3(R);
    }

    distance (other: Nucleotide, PBC=true, box?: THREE.Vector3): THREE.Vector3 {
        if (PBC && box === undefined) {
            Logger.die ("distance between nucleotides: if PBC is true, box must be provided");
        }
        let dr = other.pos_back.clone().sub(this.pos_back);
        if (PBC) {
            dr.sub(box.clone().multiply(dr.clone().divide(box).round()));
        }
        return dr;
    }

    get_base(): string {
        /*
        Returns a number containing base id
        >>> v1 = new THREE.Vector3(0.,0.,0.);
        >>> v2 = new THREE.Vector3(1.,0.,0.);
        >>> v3 = new THREE.Vector3(0.,0.,1.);
        >>> new Nucleotide(v1, v2, v3, 'A').get_base()
        'A'
        >>> new Nucleotide(v1, v2, v3, "C").get_base()
        'C'
        >>> new Nucleotide(v1, v2, v3, "G").get_base()
        'G'
        >>> new Nucleotide(v1, v2, v3, "T").get_base()
        'T'
        >>> new Nucleotide(v1, v2, v3, 1).get_base()
        'G'
        >>> new Nucleotide(v1, v2, v3, 103).get_base()
        '103'
        >>> new Nucleotide(v1, v2, v3, -97).get_base()
        '-97'
        */
        
        if ([0, 1, 2, 3].includes(this._base)) {
            return number_to_base[this._base];
        } else {
            return this._base.toString();
        }
    }

    _get_lorenzo_output() {
        return [this.cm_pos, this._a1, this._a3, this._v, this._L].map(v=>v.toArray().join(' ')).join(' ');
    }
}

class Strand {
    /*
    Strand composed of Nucleotides
    Strands can be contained in System
    */

    index: number;
    _first: number;
    _last: number;
    _nucleotides: Nucleotide[];
    _sequence: number[];
    _circular: boolean;  //  bool on circular DNA

    constructor() {
        this.index = Strand.index;
        Strand.index++;
        this._first = -1;
        this._last = -1;
        this._nucleotides = [];
        this._sequence = [];
        this._circular = false;  //  bool on circular DNA
    }
    static index = 0;

    get N() {
        return this._nucleotides.length;
    }

    get sequence() {
        return this._sequence;
    }

    _prepare(si: number, ni: number) {
        this.index = si;
        this._first = ni;

        let n=0;
        while (n<this.N) {
            this._nucleotides[n].index = ni + n;
            n++;
        }

        this._last = ni + n;
        return ni + n + 1;
    }

    copy() {
        let copy = new Strand();
        for (const n of this._nucleotides) {
            copy.add_nucleotide(n.copy())
        }
        return copy;
    }

    get cm_pos() {
        let total = new THREE.Vector3();
        this._nucleotides.forEach(n=>{
            total.add(n.cm_pos);
        })
        return total.divideScalar(this.N);
    }

    set cm_pos(new_pos) {
        const diff = new_pos.clone().sub(this.cm_pos);
        this._nucleotides.forEach(n=>n.translate(diff));
    }

    translate(amount) {
        let new_pos = this.cm_pos.clone().add(amount);
        this.cm_pos = new_pos;
    }

    rotate(R: THREE.Matrix3, origin?: THREE.Vector3) {
        if (origin === undefined) { 
            origin = this.cm_pos;
        }

        for (const n of this._nucleotides) { 
            n.rotate(R, origin);
        }
    }

    append(other: Strand) {
        let dr = this._nucleotides.slice(-1)[0].distance(other._nucleotides[0], false);
        if (Math.sqrt(dr.dot(dr)) > (0.7525 + 0.25)) {
            Logger.log("WARNING: Strand.push(): strands seem too far apart. Assuming you know what you are doing.")
        }

        let ret = new Strand();

        for (const n of this._nucleotides) {
            ret.add_nucleotide(n);
        }

        for (const n of other._nucleotides) {
            ret.add_nucleotide(n);
        }

        return ret;
    }

    get_slice(start=0, end?: number) {
        if (end === undefined) { 
            end = this.N;
        }

        if (end > this.N) {
            throw `The given end parameter is larger than the number of nucleotides of the strand (${end} > ${this.N})`;
        }

        let ret = new Strand();
        for (let i=start; i<end; i++) {
            ret.add_nucleotide(this._nucleotides[i].copy());
        }
        return ret;
    }

    set sequence(seq: string | number[]) {
        if (typeof seq === "string") {
            seq = Array.from(seq).map(c=>base_to_number[c]);
        }
        if (seq.length != this._nucleotides.length) {
            Logger.log ("Cannot change sequence: lengths don't match", Logger.WARNING)
            return;
        }
        for(let i=0; i<this._nucleotides.length; i++) {
            this._nucleotides[i]._base = seq[i];
        }
        this._sequence = seq as number[];
    }

    bring_in_box_nucleotides(box) {
        let diff = this.cm_pos.divide(box).round().multiply(box);
        for (const n of this._nucleotides) {
            n.cm_pos_box = n.cm_pos.clone().sub(diff);
        }
    }

    add_nucleotide(n: Nucleotide) {
        if (this._nucleotides.length === 0) {
            this._first = n.index;
        }
        n.strand = this.index;
        this._nucleotides.push(n);
        this._last = n.index;
        this._sequence.push(n._base);
    }

    _get_lorenzo_output() {
        let conf = this._nucleotides.map(n=>n._get_lorenzo_output()).join('\n') + '\n';

        let top = "";
        for (const n of this._nucleotides) {
            let n3: number;
            let n5: number;
            if (this._circular) {
                if (n.index === this._first) {
                    n3 = this._last;
                } else {
                    n3 = n.index - 1;
                }
                if (n.index === this._last) {
                    n5 = this._first;
                } else {
                    n5 = n.index + 1;
                }
            } else {
                if (n.index === this._first) {
                    n3 = -1;
                } else {
                    n3 = n.index - 1;
                }
                if (n.index === this._last) {
                    n5 = -1;
                } else {
                    n5 = n.index + 1;
                }
            }
            top += `${this.index+1} ${n.get_base()} ${n3} ${n5}\n`;
        }

        return [conf, top];
    }

    get_lammps_N_of_bonds_strand() {
        let N_bonds = 0;
        for (const n of this._nucleotides) {
            if (n.index != this._last) {
                N_bonds++;
            } else if (this._circular) {
                N_bonds++;
            }
        }

        return N_bonds;
    }

    get_lammps_bonds(): string[] {
        let top = [];
        for (const n of this._nucleotides) {
            if (n.index != this._last) {
                top.push(`${n.index+1}  ${n.index+2}`);
            } else if (this._circular) {
                top.push(`${n.index+1}  ${this._first + 1}`);
            }
        }
        return top;
    }

    make_circular(check_join_len=false) {
        if (check_join_len) {
            let dr = this._nucleotides.slice(-1)[0].distance(this._nucleotides[0], false);
            if (Math.sqrt(dr.dot(dr)) > (0.7525 + 0.25)) {
                Logger.log("Strand.make_circular(): ends of the strand seem too far apart. Assuming you know what you are doing.", Logger.WARNING);
            }
        }
        this._circular = true;
    }

    make_noncircular() {
        this._circular = false;
    }
        
    is_circular() {
        return this._circular;
    }

    cut_in_two(copy=true) { //  cuts a strand into two strands in the middle
        let fragment_one = new Strand();
        let fragment_two = new Strand();
        let counter = 0;
        for (const n of this._nucleotides) {
            if (counter < (this._nucleotides.length/2)) {
                fragment_one.add_nucleotide(copy? n.copy() : n);
            } else {
                fragment_two.add_nucleotide(copy? n.copy() : n);
            }
            counter++;
        }
        return [fragment_one, fragment_two];
    }
}


class System {
    /*
    Object representing an oxDNA system
    Contains strands

    Arguments:
    box -- the box size of the system
        Ex: box = [50, 50, 50];

    time --- Time of the system

    E_pot --- Potential energy

    E_kin --- Kinetic energy

    */

    _time: number;
    _ready: boolean;
    _box: THREE.Vector3;
    _N: number
    _N_strands: number;
    _strands: Strand[];
    _nucleotide_to_strand: number[];
    _N_cells: THREE.Vector3;
    _cellsides: THREE.Vector3;
    E_pot: number;
    E_kin: number;
    E_tot: number;

    constructor(box: THREE.Vector3, time=0, E_pot=0, E_kin=0) {
        this._time = time;
        this._ready = false;
        this._box = box;
        this._N = 0;
        this._N_strands = 0;
        this._strands = [];
        this._nucleotide_to_strand = [];
        this.E_pot = E_pot;
        this.E_kin = E_kin;
        this.E_tot = E_pot + E_kin;
        
        Nucleotide.index = 0;
        Strand.index = 0;
    }

    get sequences() {
        return this._strands.map(x=>x._sequence);
    }

    get N() {
        return this._N;
    }

    get N_strands() {
        return this._N_strands;
    }

    _prepare() {
        let nind = 0;
        for (let sind=0; sind<this._N_strands; sind++) {
            nind = this._strands[sind]._prepare(sind, nind);
        }

        for (const s of this._strands) {
            s.bring_in_box_nucleotides(this._box);
        }
    }

    copy () {
        let copy = new System(this._box);
        for (const s of this._strands) {
            copy.add_strand(s.copy ());
        }
        return copy;
    }

    join(other: System, box: THREE.Vector3) {
        if (box === undefined) {
            box = new THREE.Vector3(0., 0., 0.);
            for(let i=0; i<3; i++) {
                if (other._box[i] > this._box[i]) {
                    box[i] = other._box[i];
                } else {
                    box[i] = this._box[i];
                }
            }
        }

        let ret = new System(box);
        for (const s of this._strands) {
            ret.add_strand(s.copy());
        }
        for (const s of other._strands) {
            ret.add_strand(s.copy());
        }

        return ret;
    }

    add_strand(s: Strand) {
        // Add a Strand to the System
        this._strands.push(s);
        this._N += s.N;
        this._N_strands++;
        return true;
    }

    add_strands(ss: Strand | Strand[]) {
        if (Array.isArray(ss)) {
            let added = [];
            for (const s of ss) {
                if (this.add_strand(s)) {
                    added.push(s)
                }
            }
            if (added.length === ss.length) {
                return true
            } else {
                for (const s of added) {
                    Nucleotide.index -= s.N;
                    Strand.index--;
                    this._strands.pop();
                    this._N -= s.N;
                    this._N_strands--;
                }
                return false;
            }
        } else if (!this.add_strand(ss)) { 
            return false;
        }
        return true;
    }

    rotate (amount, origin) {
        for (const s of this._strands) {
            s.rotate (amount, origin)
        }
    }

    translate (amount) {
        for (const s of this._strands) {
            s.translate (amount);
        }
    }

    print_lorenzo_output(): [string, string] {
        let conf = `t = ${this._time}\nb = ${this._box.x} ${this._box.y} ${this._box.z}\nE = ${this.E_tot} ${this.E_pot} ${this.E_kin}\n`;

        let visible_strands = 0;
        let visible_nucleotides = 0;
        for (const s of this._strands) {
            visible_strands++;
            visible_nucleotides += s.N;
        }

        let topology = `${visible_nucleotides} ${visible_strands}\n`;
        for (const s of this._strands) {
            let [sc, st] = s._get_lorenzo_output();
            topology += st;
            conf += sc;
        }

        return [topology, conf];
    }

    print_oxview_output() {
        let out = {
            'box': this._box.round().toArray(),
            'systems': [{'id':0, 'strands': []}]
        }

        for (const s of this._strands) {
            let strand = {
                'id': s.index, 'end3': s._nucleotides[0].index, 'end5': s._nucleotides.slice(-1)[0].index,
                'class': 'NucleicAcidStrand', 'monomers': []
            };
            for(let i=0; i<s.N; i++) {
                let n = s._nucleotides[i];
                let n5: number;
                let n3: number;
                if (s._circular) {
                    if (i === 0) {
                        n3 = s._nucleotides.slice(-1)[0].index;
                    } else {
                        n3 = s._nucleotides[i-1].index;
                    }
                    if (i === s._nucleotides.length-1) {
                        n5 = s._nucleotides[0].index;
                    } else {
                        n5 = s._nucleotides[i+1].index;
                    }
                } else {
                    if (i === 0) {
                        n3 = -1;
                    } else {
                        n3 = s._nucleotides[i-1].index;
                    }
                    if (i === s._nucleotides.length-1) {
                        n5 = -1;
                    } else {
                        n5 = s._nucleotides[i+1].index;
                    }
                }
                let nucleotide = {
                    'id': n.index,
                    'type': n.get_base(),
                    'class': 'DNA',
                    'p': n.cm_pos.toArray(),
                    'a1': n._a1.toArray(),
                    'a3': n._a3.toArray()
                }

                if (n3 >= 0) nucleotide['n3'] = n3;
                if (n5 >= 0) nucleotide['n5'] = n5;
                if (n.pair !== undefined) nucleotide['bp'] = n.pair.index;
                if (n.cluster !== undefined) nucleotide['cluster'] = n.cluster;
                if (n.color !== undefined) nucleotide['color'] = n.color;

                strand['monomers'].push(nucleotide);
            }
            out['systems'][0]['strands'].push(strand);
        }

        return JSON.stringify(out);
    }

    get _nucleotides () {
        return [].concat(...this._strands.map(s=>s._nucleotides));
    }

    map_nucleotides_to_strands() {
        //  this function creates nucl_id -> strand_id array
        for(let i=0; i<this._strands.length; i++) {
            for (let j=0; j<this._strands[i].N; j++) {
                this._nucleotide_to_strand.push(i);
            }
        }
    }

    print_dot_bracket_output() {
        //  assumes each nucleotide has at most 1 hydrogen bond, requires interactions already to be filled for nucleotide objects
        let nupack_string = "";
        for (let n1=0; n1<this.N; n1++) {
            let interactions = this._nucleotides[n1].interactions;
            if (interactions.length > 1) {
                Logger.log ("more than 1 HB for a nucleotide", Logger.WARNING);
            }
            if (interactions.length === 0) {
                nupack_string += "."
            } else if (interactions[0] > n1) {
                nupack_string += "("
            } else if (interactions[0] < n1) {
                nupack_string += ")"
            } else {
                Logger.log("unexpected interaction detected while building nupack string", Logger.CRITICAL)
            }
        }
        return nupack_string;
    }
}

export {System, Strand, Nucleotide, Logger,
    base_to_number, FLT_EPSILON,
    POS_BACK, POS_STACK, POS_BASE, CM_CENTER_DS,
    FENE_R0_OXDNA, FENE_EPS,LENGTH_FACT, BASE_BASE,
    INT_HYDR, INT_STACK,INT_CROSS_STACK, INT_COAX_STACK, INT_FENE, INT_EXC_BONDED, INT_EXC_NONBONDED,
    H_CUTOFF, MM_GROOVING,
};