/*
Created on Nov 11, 2018

@author: lorenzo
*/


import * as base from "./base"
import * as utils from "./utils"
import * as THREE from 'three'


class StrandGenerator {

    generate(bp, sequence=undefined, start_pos=new THREE.Vector3(0, 0, 0), direction=new THREE.Vector3(0, 0, 1), perp=undefined, rot=0., double=true, circular=false, DELTA_LK=0, BP_PER_TURN=10.34, ds_start=undefined, ds_end=undefined, force_helicity=false) {
        /*
        Generate a strand of DNA.
            - linear, circular (circular)
            - ssDNA, dsDNA (double)
            - Combination of ss/dsDNA (ds_start, ds_end)
            Note: Relevent argument(s) in parentheses.

        Arguments:
        bp --- Integer number of bp/nt (required)
        sequence --- Array of integers or string. Should be same length as bp (default undefined)
            Default (undefined) generates a random sequence.
            Ex: [0,1,2,3,0]
            Ex: "AGCTA"
            See dictionary base.base_to_number for int/char conversion {0:'A'}
        start_pos --- Location to begin building the strand (default new THREE.Vector3(0, 0, 0))
        direction --- a3 vector, orientation of the base (default new THREE.Vector3(0, 0, 1))
        perp --- Sets a1 vector, the orientation of the backbone. (default false)
            Must be perpendicular to direction (as a1 must be perpendicular to a3)
            If perp === undefined or false, perp is set to a random orthogonal angle
        rot --- Rotation of first bp (default 0.)
        double --- Generate dsDNA (default true)
        circular --- Generate closed circular DNA (defalt false)
            Limitations...
            For ssDNA (double=false): bp >= 4
            For dsDNA (double=true) : bp >= 30
            Will throw warnings. Allowed, but use at your own risk.
        DELTA_LK --- Integer change in linking number from Lk0 (default 0)
            Only valid if circular===true
        BP_PER_TURN --- Base pairs per complete 2*pi helix turn. (default 10.34)
            Only valid if circular===true
        ds_start --- Index (from 0) to begin double stranded region (default undefined)
        ds_end --- Index (from 0) to end double stranded region (default undefined)
            Default === undefined, which is entirely dsDNA; sets ds_start = 0, ds_end=bp;
            Ex: ds_start=0, ds_end=10 will create a double stranded region on bases
                range(0,10): [0,1,2,3,4,5,6,7,8,9]
            Note: To generate a nicked circular dsDNA, manually change state with
                  {Strand}.make_noncircular()
        force_helicity --- Force generation of helical strands. Use helicity by default
            for bp > 30. Warns from 18 to 29. Will crash oxDNA below 18. (default false)

        Note: Minimuim circular duplex is 18. Shorter circular strands disobey FENE.
        For shorter strands, circular ssDNA is generated in a circle instead of having
        imposed helicity.

        Examples:
        Generate ssDNA:
            generate(bp=4,sequence=[0,1,2,3],double=false,circular=false)
        Generate circular dsDNA with +2 Linking number:
            generate(bp=45,double=true,circular=true,DELTA_LK=2)
        Generate a circular ssDNA (45nt) with ssDNA (25nt) annealed to indices 0 to 24:
            generate(bp=45,double=true,circular=true,ds_start=0,ds_end=25)
        */

        //  Loads of input checking...
        if (typeof sequence === 'string') {
            try {
                sequence = Array.from(sequence).map(c=>base.base_to_number[c]);
            }
            catch (e) {
                base.Logger.die("Key Error: sequence is invalid")
            }
        }
        if (sequence === undefined) {
            sequence = new Array(bp);
            for(let i=0; i<bp; i++) {
                sequence[i] = base.base_to_number['T'];
            }
        } else if (sequence.length != bp) {
            let n = bp - sequence.length;
            base.Logger.log(`sequence is too short, adding ${n} thymine bases`, base.Logger.WARNING);
            while(n--) {
                sequence.push(
                    base.base_to_number['T']
                )
            }
        }

        if (circular && bp < 30) {
            //  30 is about the cut off for circular dsDNA. Anything shorter will probably clash.
            //  oxDNA can relax down to 18.
            //  4 is about the cut off for circular ssDNA. Use dsDNA cutoff for saftey.
            base.Logger.log("sequence is too short! Proceed at your own risk", base.Logger.WARNING)
        }

        let option_use_helicity = true;
        if (circular && bp < 30 && !double) {
            base.Logger.log("sequence is too short! Generating ssDNA without imposed helicity", base.Logger.WARNING)
            //  Do not impose helcity to generate shorter circular ssDNA
            if (!force_helicity) {
                option_use_helicity = false;
            }
        }

        if (ds_start === undefined) {
            ds_start = 0;
        }
        if (ds_end === undefined) {
            ds_end = bp;
        }
        if (ds_start > ds_end) {
            base.Logger.die("ds_end > ds_start")
        }
        if (ds_end > bp) {
            base.Logger.die("ds_end > bp")
        }

        //  we need to find a vector orthogonal to direction
        let dir_norm = Math.sqrt(direction.dot(direction));
        if (dir_norm < 1e-10) {
            base.Logger.log("direction must be a valid vector, defaulting to (0, 0, 1)", base.Logger.WARNING)
            direction = new THREE.Vector3(0, 0, 1);
        } else {
            direction.divideScalar(dir_norm);
        }
        let v1: THREE.Vector3;
        if (!perp) {
            v1 = new THREE.Vector3(
                Math.random(), Math.random(), Math.random()
            );
            v1.sub(direction.clone().multiplyScalar(direction.clone().dot(v1)));
            v1.normalize();
        } else {
            v1 = perp;
        }

        //  Setup initial parameters
        let ns1 = new base.Strand();
        //  and we need to generate a rotational matrix

        let a1 = v1.clone().applyAxisAngle(direction, rot);
        let rb = start_pos.clone()
        let a3 = direction;

        let torus_perp: THREE.Vector3, angle: number, radius: number;
        //  Circular strands require a continuious deformation of the ideal helical pitch
        if (circular) {
            //  Unit vector orthogonal to plane of torus
            //  Note: Plane of torus defined by v1,direction
            torus_perp = v1.clone().cross(direction);
            //  Angle between base pairs along torus
            angle = 2. * Math.PI / bp;
            //  Radius of torus
            radius = base.FENE_R0_OXDNA / Math.sqrt(2. * (1. - Math.cos(angle)));
        }

        if (circular && option_use_helicity) {
            //  Draw backbone in a helical spiral around a torus
            //  Draw bases pointing to center of torus
            for (let i=0; i<bp; i++) {
                //  Torus plane defined by direction and v1
                let v_torus = v1.clone().multiplyScalar(base.BASE_BASE * Math.cos(i * angle)).add(direction.clone().multiplyScalar(base.BASE_BASE * Math.sin(i * angle)));
                rb.add(v_torus);

                //  a3 is tangent to the torus
                a3 = v_torus.clone().normalize();

                //  a1 is orthogonal to a3 and the torus normal
                a1 = a3.clone().cross(torus_perp);

                //  Apply the rotation matrix
                a1.applyAxisAngle(a3, i * (Math.round(Math.floor(bp / BP_PER_TURN)) + DELTA_LK) / bp * 2* Math.PI)
                ns1.add_nucleotide(new base.Nucleotide(rb.clone().sub(a1.clone().multiplyScalar(base.CM_CENTER_DS)), a1, a3, sequence[i]));
            }
            ns1.make_circular(true);
        } else if (circular && !option_use_helicity) {
            for(let i=0; i<bp; i++) {
                rb = new THREE.Vector3(
                    Math.cos (i * angle) * radius + 0.34 * Math.cos(i * angle),
                    Math.sin (i * angle) * radius + 0.34 * Math.sin(i * angle),
                    0,
                );
                a1 = new THREE.Vector3(
                    Math.cos (i * angle),
                    Math.sin (i * angle),
                    0
                );
                ns1.add_nucleotide(new base.Nucleotide(rb, a1, new THREE.Vector3(0, 0, 1), sequence[i]));
            }
            ns1.make_circular(true);
        } else {
            //  Add nt in canonical double helix
            for(let i=0; i<bp; i++) {
                ns1.add_nucleotide(new base.Nucleotide(rb.clone().sub(a1.clone().multiplyScalar(base.CM_CENTER_DS)), a1, a3, sequence[i]))
                if (i != bp - 1) {
                    a1.applyAxisAngle(direction, 35.9 * Math.PI/180);
                    rb.add(a3.multiplyScalar(base.BASE_BASE));
                }
            }
        }

        //  Fill in complement strand
        if (double) {
            let ns2 = new base.Strand();
            for (let i=ds_end-1; i<=ds_start; i--) {
                //  Note that the complement strand is built in reverse order
                let nt = ns1._nucleotides[i];
                a1 = nt._a1.clone().negate();
                a3 = nt._a3.clone().negate();
                let nt2_cm_pos = a1.clone().multiplyScalar(-(base.FENE_EPS + 2 * base.POS_BACK)).add(nt.cm_pos);
                ns2.add_nucleotide(new base.Nucleotide(nt2_cm_pos, a1, a3, 3 - sequence[i]));
            }
            if (ds_start === 0 && ds_end === bp && circular) {
                ns2.make_circular(true);
            }
            return [ns1, ns2]
        } else {
            return ns1
        }
    }

    generate_or_sq(bp: number, start_pos=new THREE.Vector3(0., 0., 0.), direction=new THREE.Vector3(0., 0., 1.), perp=undefined, double=true, rot=0., angle: number | number[] = Math.PI / 180 * 33.75, length_change=[], region_begin=[], region_end=[]) {
        if (length_change && region_begin.length != region_end.length) {
            if ((region_end.length + 1) === region_begin.length) {
                base.Logger.log(`the lengths of begin ${region_begin.length} and end ${region_end.length} arrays are mismatched; I will try to proceed by using the number of basepairs as the last element of the end array`, base.Logger.WARNING)
                region_end.push(bp + 1)
            } else {
                base.Logger.die(`the lengths of begin ${region_begin.length} and end ${region_end.length} arrays are unrecoverably mismatched`);
            }
        }

        //  angle should be an array, with a length 1 less than the # of base pairs
        if (!Array.isArray(angle)) {
            let tmp = angle;
            angle = [];
            for(let i=0; i<bp; i++) {
                angle.push(tmp);
            }
        } else if (angle.length != bp - 1) {
            base.Logger.log("generate_or_sq: incorrect angle array length, should be 1 less than number of base pairs", base.Logger.CRITICAL);
        }

        //  we need to find a vector orthogonal to direction
        let dir_norm = Math.sqrt(direction.dot(direction));
        if (dir_norm < 1e-10) {
            base.Logger.log("direction must be a valid vector, defaulting to (0, 0, 1)", base.Logger.WARNING)
            direction = new THREE.Vector3(0, 0, 1);
        } else { 
            direction.divideScalar(dir_norm);
        }

        let v1: THREE.Vector3;
        if (perp === undefined) {
            v1 = new THREE.Vector3(Math.random(), Math.random(), Math.random());
            v1.sub(direction.clone().multiplyScalar(direction.clone().dot(v1)))
            v1.normalize();
        } else {
            v1 = perp.clone();
        }

        let ns1 = new base.Strand();
        let a1: THREE.Vector3 = v1.clone().applyAxisAngle(direction, rot);
        let rb: THREE.Vector3 = start_pos.clone();
        let a3: THREE.Vector3 = direction.clone();
        let Rs = [];

        for(let i=0; i<bp; i++) {
            ns1.add_nucleotide(new base.Nucleotide(
                rb.clone().sub(a1.clone().multiplyScalar(base.CM_CENTER_DS)),
                a1.clone(), a3.clone(), undefined)
            );
            if (i != bp - 1) {
                let R = new THREE.Quaternion().setFromAxisAngle(direction, angle[i]);
                Rs.push(R);
                a1.applyQuaternion(R).normalize();
                rb.add(a3.clone().normalize().multiplyScalar(base.BASE_BASE));
                if (length_change) {
                    for (let j=0; j<length_change.length; j++) {
                        if (i >= region_begin[j] && i < region_end[j]) {
                            if (length_change[j]) {
                                rb.add(a3.clone().multiplyScalar(base.BASE_BASE * (-(length_change[j]) / (region_end[j] - region_begin[j]))));
                            }
                        }
                    }
                }
            }
        }

        if (double) {
            a1.negate().normalize();
            a3 = direction.clone().negate().normalize();
            let ns2 = new base.Strand();

            for(let i=0; i<bp; i++) {
                //  create new nucleotide and save basepair info on both sides
                let paired_nuc = ns1._nucleotides[bp-i-1];
                let new_nuc = new base.Nucleotide(
                    rb.clone().sub(a1.clone().multiplyScalar(base.CM_CENTER_DS)),
                    a1.clone(), a3.clone(), undefined, undefined, undefined,
                    undefined, paired_nuc
                );
                paired_nuc.pair = new_nuc;
                ns2.add_nucleotide(new_nuc)
                if (i != bp - 1) {
                    //  we loop over the THREE.Quaternions in the reverse order, and use the inverse of rotation
                    a1.applyQuaternion(Rs.pop().conjugate()).normalize();
                    rb.add(a3.clone().multiplyScalar(base.BASE_BASE))
                    if (length_change) {
                        for (let j=0; j<length_change.length; j++) {
                            if (bp - 2 - i >= region_begin[j] && bp - 2 - i < region_end[j]) {
                                if (length_change[j]) {
                                    rb.add(a3.clone().multiplyScalar(base.BASE_BASE * (-(length_change[j] / (region_end[j] - region_begin[j])))))
                                }
                            }
                        }
                    }
                }
            }

            return [ns1, ns2];
        } else {
            return ns1;
        }
    }

    generate_double_offset(seqA: string | number[], seqB: string | number[], offset, start_pos=new THREE.Vector3(0, 0, 0), direction=new THREE.Vector3(0, 0, 1), perp=undefined, rot=0): base.Strand[] {
        let seqa: number[];
        if (typeof seqA === "string") {
            seqa = [...seqA].map(x=>base.base_to_number[x]);
        } else {
            seqa = seqA;
        }
        let seqb: number[];
        if (typeof seqB === "string") {
            seqb = [...seqB].map(x=>base.base_to_number[x]);
        } else {
            seqb = seqB;
        }

        let bp = Math.max(seqa.length, seqb.length + offset);

        let strands = this.generate(bp, undefined, start_pos, direction, false, undefined, true);
        let s1 = strands[0];
        let s2 = strands[1];

        s1 = s1.get_slice (0, seqa.length);

        if (seqb.length + offset > seqa.length) {
            s2 = s2.get_slice (0, seqb.length)  //  starts from opposite end;
        } else {
            s2 = s2.get_slice (bp - offset - seqb.length, seqb.length);
        }
        s1.sequence = seqa;
        s2.sequence = seqb;

        return [s1, s2]
    }

    generate_rw (sequence, start_pos=new THREE.Vector3(0., 0., 0.)) {
        /*
        Generate ssDNA as a random walk (high-energy configurations are possible):
            generate(bp=45,double=false,circular=false,random_walk=true)
        */
        //  random walk generator
        base.Logger.log("Generating strand as a random walk. Remember to equilibrate the configuration with MC", base.Logger.WARNING)
        let d = new THREE.Vector3(0.7525, 0., 0.);
        let pos = start_pos;
        let rw = [];
        rw.push(pos)
        for (let i=1; i<sequence.length; i++) {
            let overlap = true;
            let trypos: THREE.Vector3;
            while (overlap) {
                overlap = false;
                let R = utils.get_random_rotation_matrix();
                let dd = d.clone().applyMatrix3(R);
                trypos = pos.clone().add(d.clone().applyMatrix3(R));
                overlap = false;
                for (const r of rw) {
                    dd = trypos.clone().sub(r);
                    if (dd.dot(dd) < 0.40 * 0.40) {
                        overlap = true;
                    }
                }
            }
            pos = trypos;
            rw.push (pos)
        }

        //  we get the a1 vectors in a smart way
        let a1s = [];
        d = rw[1].clone().sub(rw[0]);
        a1s.push(d.clone().divideScalar(Math.sqrt (d.dot(d))))

        for (let i=1; i<rw.length-1; i++) {
            d = (rw[i + 1].clone().add(rw[i - 1])).multiplyScalar(0.5);
            d = rw[i].clone().sub(d);
            a1s.push(d.clone().divideScalar(Math.sqrt(d.dot(d))));
        }

        d = rw[rw.length - 1].clone().sub(rw[rw.length - 2]);
        a1s.push (d.divideScalar(Math.sqrt(d.dot(d))));

        let s = new base.Strand();
        for (let i=0; i<rw.length; i++) {
            let r = rw[i];
            let [a1, _, a3] = utils.get_orthonormalized_base (a1s[i], utils.get_random_vector(), utils.get_random_vector()) ;
            //  we use abs since POS_BACK is negative
            let cm = r + a1s[i] * Math.abs(base.POS_BACK);
            s.add_nucleotide(new base.Nucleotide (cm, a1, a3, sequence[i]));
        }

        return s;
    }
}

// Used since array equality doesn't work in javascript
class PairMap {
    map: Map<string,any>;
    constructor() {
        this.map = new Map();
    }
    set(key: [number, number] , val: any) {
        this.map.set(key.toString(), val);
    }
    get(key: [number, number]) {
        return this.map.get(key.toString());
    }

    has(key: [number, number]) {
        return this.map.has(key.toString());
    }

    get size() {
        return this.map.size;
    }

    *keys() {
        for (let k of this.map.keys()) {
            yield k.split(',').map(v=>parseInt(v));
        }
    }
    *entries() {
        for (let [k,v] of this.map.entries()) {
            let keyPair = k.split(',').map(v=>parseInt(v)) as [number, number];
            yield [keyPair, v];
        }
    }
    values() {
        return this.map.values();
    }
}

class vhelix_vbase_to_nucleotide extends PairMap {
    _scaf: PairMap;
    _stap: PairMap;
    nuc_count: number;
    strand_count: number;
    //  at the moment squares with skips in have entries in the dicts but with the nucleotide list empty (rather than having no entry) - I'm not sure whether or not this is desirable. It's probably ok
    constructor() {
        super();
        this._scaf = new PairMap();
        this._stap = new PairMap();
        this.nuc_count = 0 //  record the nucleotide count, updated only after a whole strand is added;
        this.strand_count = 0;
    }

    add_scaf(vh, vb, strand, nuc) {
        this._scaf.set([vh, vb], [strand, nuc]);
    }

    add_stap(vh, vb, strand, nuc) {
        this._stap.set([vh, vb], [strand, nuc]);
    }

    //  these methods use a reference vhvb2n object to make the final vhvb2n object
    add_scaf_strand(add_strand, reference, continue_join = false) {;
        let count = 0;
        const size = this._scaf.size;
        for (const [[vh, vb], [strand_ind, nuc]] of reference._scaf.entries()) {
            if (strand_ind === add_strand) {
                this.add_scaf(vh, vb, this.strand_count, nuc.map(x=>x + this.nuc_count));
                count += nuc.length;
            }
        }
        this.nuc_count += count;
        if (this._scaf.size === size) {
            return 1;
        } else {
            if (!continue_join) {
                this.strand_count++;
            }
            return 0;
        }
    }

    add_stap_strand(add_strand, reference, continue_join = false) {;
        let count = 0;
        const size = this._stap.size;
        for (const [[vh, vb], [strand_ind, nuc]] of reference._stap.entries()) {
            if (strand_ind === add_strand) {
                this.add_stap(vh, vb, this.strand_count, nuc.map(x=>x + this.nuc_count));
                count += nuc.length;
            }
        }
        this.nuc_count += count;
        if (this._stap.size === size) {
            return 1;
        } else {
            if (!continue_join) {
                this.strand_count++;
            }
            return 0;
        }
    }

    add_strand(add_strand, reference, continue_join = false) {;
        return (this.add_scaf_strand(add_strand, reference, continue_join) && this.add_stap_strand(add_strand, reference, continue_join));
    }
}

export {StrandGenerator, vhelix_vbase_to_nucleotide, PairMap}
