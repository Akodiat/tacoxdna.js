import * as THREE from 'three'
import * as base from "./libs/base";
import * as cu from "./libs/cadnano_utils";
import * as utils from "./libs/utils";

const DIST_HEXAGONAL = 2.55  //  distance between centres of virtual helices (hexagonal array);
const DIST_SQUARE = 2.60  //  distance between centres of virtual helices (square array);
const BOX_FACTOR = 2  //  factor by which to expand the box (linear dimension);


class vh_nodes {
    begin: number[];
    end: number[];
    constructor() {
        this.begin = [];
        this.end = [];
    }

    toString() {
        return `${this.begin} ${this.end}`;
    }

    add_begin(begin_index: number) {
        if (!this.begin.includes(begin_index)) {
            this.begin.push(begin_index);
        }
    }

    add_end(end_index: number) {
        if (!this.end.includes(end_index)) {
            this.end.push(end_index)
        }
    }
}

function vhelix_rotation_origami_sq(direction: THREE.Vector3, perp: THREE.Vector3) {
    return perp.clone().applyAxisAngle(direction, Math.PI * 15. / 180);
}


function vhelix_rotation_origami_he(direction: THREE.Vector3, perp: THREE.Vector3) {
    return perp.clone().applyAxisAngle(direction, Math.PI * 160. / 180);
}

function insert_loop_skip(start_pos, direction, perp, rot, helix_angles, vhelix: vhelix, nodes: vh_nodes) {
    //  return a double strand which is a copy of the double strand in the first argument, but with skips and loops

    //  strand is generated right to left i.e. opposite direction to even vhelix
    let length_change = [];
    let length_change_total = 0;
    let new_nodes = new vh_nodes();
    let new_angle = [];
    let helix_angles_new: number[] = helix_angles.slice(); // copy

    let reverse_nodes: vh_nodes;
    if (vhelix.num % 2 === 1) {
        reverse_nodes = new vh_nodes();
        reverse_nodes.begin = nodes.begin.slice().reverse();
        reverse_nodes.end = nodes.end.slice().reverse();
    }
    for(let i=0; i<nodes.begin.length; i++) {
        //  ltr: left to right; looking at the strand left to right (low to high square index), the beginning/end of the effective strand is here (before skips/loops)
        //  gs: generated strand; the index of the nucleotide (BEFORE skips/loops are applied) on the generated strand corresponding to the beginning/end of the effective strand
        let begin_ltr, begin_gs, end_gs, end_ltr;
        if (vhelix.num % 2 === 0) {
            begin_ltr = nodes.begin[i];
            end_ltr = nodes.end[i];
            begin_gs = nodes.begin[i];
            end_gs = nodes.end[i];
        } else {
            begin_ltr = reverse_nodes.end[i];
            end_ltr = reverse_nodes.begin[i];
            begin_gs = vhelix.len - reverse_nodes.begin[i] - 1;
            end_gs = vhelix.len - reverse_nodes.end[i] - 1;
        }

        //  check for zero length effective strand
        if (end_gs - begin_gs != 0) {
            //  get length change for this effective strand
            length_change.push(0)
            for (const j of vhelix.skip.slice(begin_ltr,end_ltr + 1)) {
                length_change[i] -= j;
            }
            for (const j of vhelix.loop.slice(begin_ltr, end_ltr + 1)) {
                length_change[i] += j;
            }
            //  get new pitch angles for this effective strand
            new_angle.push(utils.sum(...helix_angles.slice(begin_gs, end_gs)) / (end_gs - begin_gs + length_change[i]))
            for (let j=begin_gs; j<end_gs; j++) {
                helix_angles_new[j] = new_angle[i];
            }
            //  adjust beginning/end indices according to length change
            begin_gs += length_change_total;
            end_gs += length_change_total + length_change[i];
            new_nodes.add_begin(begin_gs);
            new_nodes.add_end(end_gs);  //  begin_gs > end_gs.....
        } else {
            length_change.push(0);
            new_angle.push(utils.sum(...helix_angles) / helix_angles.length);  //  append an average angle
            //  adjust beginning/end indices according to length change
            begin_gs += length_change_total;
            end_gs += length_change_total + length_change[i];
            new_nodes.add_begin(begin_gs);
            new_nodes.add_end(end_gs);
        }
        length_change_total += length_change[i];
    }

    //  adjust the new helix angle array according to skips/loops
    let deleted = 0;
    let inserted = 0;
    let deleted_this_iteration = 0;
    let inserted_this_iteration = 0;
    for(let i=0; i<nodes.begin.length; i++) {
        deleted += deleted_this_iteration
        inserted += inserted_this_iteration
        deleted_this_iteration = 0;
        inserted_this_iteration = 0;
        let begin_ltr, end_ltr, begin_gs, end_gs;
        if (vhelix.num % 2 === 0) {
            begin_ltr = nodes.begin[i];
            end_ltr = nodes.end[i];
            begin_gs = nodes.begin[i];
            end_gs = nodes.end[i];
        } else {
            begin_ltr = reverse_nodes.end[i];
            end_ltr = reverse_nodes.begin[i];
            begin_gs = vhelix.len - reverse_nodes.begin[i] - 1;
            end_gs = vhelix.len - reverse_nodes.end[i] - 1;
        }
        for (const j of (vhelix.skip.slice(begin_ltr, end_ltr + 1))) {
            if (j === 1) {
                helix_angles_new.splice(begin_gs - deleted + inserted, 1);
                deleted_this_iteration++;
            }
        }
        for (const j of vhelix.loop.slice(begin_ltr,end_ltr + 1)) {
            for (let jj=0; jj<j; jj++) {
                helix_angles_new.splice(begin_gs - deleted + inserted, 0, new_angle[i]);
                inserted_this_iteration++;
            }
        }
    }

    let g = new cu.StrandGenerator();
    let new_strands = g.generate_or_sq(helix_angles_new.length+1, start_pos,direction,perp, true, rot, helix_angles_new, length_change, new_nodes.begin, new_nodes.end);

    return new_strands;
}


function add_slice(current_system: base.System, vhelix: vhelix, begin: number, end: number, nodes, strands, pos, direction, perp, rot, helix_angles, strand_type) {
    //  add a slice of the virtual helix to the slice system, taking into account skips and loops
    let length_change_begin = 0;
    let length_change_end = 0;

    let begin_slice, end_slice;
    
    if ((vhelix.num % 2 + strand_type) % 2 === 0) {  //  strand and even num or staple and odd num
        for (const i of vhelix.skip.slice(0, begin)) {
            length_change_begin -= i;
        }
        for (const i of vhelix.skip.slice(0, end + 1)) {
            length_change_end -= i;
        }
        for (const i of vhelix.loop.slice(0,begin)) {
            length_change_begin += i;
        }
        for (const i of vhelix.loop.slice(0, end + 1)) {
            length_change_end += i;
        }
    
        begin_slice = begin + length_change_begin;
        end_slice = end + 1 + length_change_end;

    } else {
        for (const i of vhelix.skip.slice(end)) {
            length_change_end -= i;
        }
        for (const i of vhelix.skip.slice(begin + 1)) {
            length_change_begin -= i;
        }
        for (const i of vhelix.loop.slice(end)) {
            length_change_end += i;
        }
        for (const i of vhelix.loop.slice(begin + 1)) {
            length_change_begin += i;
        }

        begin_slice = vhelix.len - begin - 1 + length_change_begin;
        end_slice = vhelix.len - end + length_change_end;
    }

    let new_strands = insert_loop_skip(pos, direction, perp, rot, helix_angles, vhelix, nodes);
    current_system.add_strand(new_strands[strand_type].get_slice(begin_slice, end_slice));
    return current_system;
}


function add_slice_nupack(vhelix: vhelix, strand_number: number, begin_helix: number, end_helix: number, index_lookup: cu.vhelix_vbase_to_nucleotide, strand_type: number) {
    let length_change = 0;
    let skips = 0;
    let loops = 0;
    if ((vhelix.num % 2 + strand_type) % 2 === 0 ) {  //  strand and even num or staple and odd num
        for (const i of vhelix.skip.slice(begin_helix, end_helix + 1)) {
            length_change -= i;
        }
        for (const i of vhelix.loop.slice(begin_helix, end_helix + 1)) {
            length_change += i;
        }
    } else {
        for (const i of vhelix.skip.slice(end_helix, begin_helix + 1)) {
            length_change -= i;
        }
        for (const i of vhelix.loop.slice(end_helix, begin_helix + 1)) {
            length_change += i;
        }
    }
    let iter_length: number;
    if ((strand_type + vhelix.num % 2) % 2 === 0) {
        iter_length = end_helix - begin_helix + 1 + length_change;
    } else {
        iter_length = begin_helix + 1 - end_helix + length_change;
    }
    let nucleotide = 0;
    while (nucleotide < iter_length) {
        let vhelix_base: number;
        if ((strand_type + vhelix.num % 2) % 2 === 0) {
            vhelix_base = nucleotide + begin_helix + skips - loops;
        } else {
            vhelix_base = begin_helix - nucleotide - skips + loops;
        }
        if (vhelix.skip[vhelix_base] != 1) {
            let add_nuc;
            if ((strand_type + vhelix.num % 2) % 2 === 0) {
                add_nuc = utils.range(vhelix.loop[vhelix_base] + 1).map(x=>(nucleotide + x));
            } else {
                add_nuc = utils.range(vhelix.loop[vhelix_base] + 1).reverse().map(x=>(nucleotide + x));
            }
            if (strand_type === 0) {
                index_lookup.set([vhelix.num, vhelix_base], [strand_number, [nucleotide]]);
            } else if (strand_type === 1) {
                index_lookup.set([strand_number, nucleotide], [vhelix.num, [vhelix_base]]);
            } else if (strand_type === 2) {
                index_lookup.add_scaf(vhelix.num, vhelix_base, strand_number, add_nuc)
            } else if (strand_type === 3) {
                index_lookup.add_stap(vhelix.num, vhelix_base, strand_number, add_nuc)
            }
            nucleotide += 1 + vhelix.loop[vhelix_base];
            loops += vhelix.loop[vhelix_base];
        } else {
            if (strand_type === 2) {
                index_lookup.add_scaf(vhelix.num, vhelix_base, strand_number, []);
            } else if (strand_type === 3) {
                index_lookup.add_stap(vhelix.num, vhelix_base, strand_number, []);
            }
            skips++;
        }
    }

    return index_lookup;
}


function build_nodes(vh) {
    //  returns a vh_nodes object which contains the beginning and end square indices of each effective strand. Effective strands
    //  on a given vhelix are a combination of information on staple and scaffold strands. Together they tell us which
    //  nucleotides need to be held constant when we alter the angles/base-base distances for a section of nucleotides along a final strand
    let nodes = new vh_nodes();
    let direction = -1;
    if (vh.num % 2 === 0) {
        direction = 1;
    }
    for(let i=0; i<vh.scaf.length; i++) {
        //  need to consider what happens when I add an index to the node list that doesn't fall within the utils.range of square indices in the vhelix
        let previd = i - 1 * direction;
        let nextid = i + 1 * direction;
        let prev, prev_stap;
        if (utils.range(vh.scaf.length).includes(previd)) {
            prev = vh.scaf[previd].type(vh, previd);
            prev_stap = vh.stap[previd].type(vh, previd);
        } else {
            prev = false;
            prev_stap = false;
        }
        let next_, next_stap;
        if (utils.range(vh.scaf.length).includes(nextid)) {
            next_ = vh.scaf[nextid].type(vh, nextid);
            next_stap = vh.stap[nextid].type(vh, nextid);
        } else {
            next_ = false;
            next_stap = false;
        }
        //  now build the effective strand vh_nodes object
        if (!(
            (prev === prev_stap && (prev === 'begin' || prev === 'end')) || 
            (next_ === next_stap && (next_ === 'begin' || next_ === 'end'))
        )) {
            if (vh.scaf[i].type(vh, i) === 'empty') {
                if (vh.stap[i].type(vh, i) === 'begin') {
                    nodes.add_end(i);
                } else if (vh.stap[i].type(vh, i) === 'end') {
                    nodes.add_begin(i);
                }
            } else if (vh.scaf[i].type(vh, i) === 'begin') {
                if (vh.stap[i].type(vh, i) === 'empty') {
                    nodes.add_begin(i);
                } else if (vh.stap[i].type(vh, i) === 'continue') {
                    nodes.add_begin(i);
                    nodes.add_end(i - 1 * direction);
                } else if (vh.stap[i].type(vh, i) === 'begin') {
                    nodes.add_begin(i + 1 * direction);
                    nodes.add_end(i - 1 * direction);
                } else if (vh.stap[i].type(vh, i) === 'end') {
                    nodes.add_begin(i);
                }
            } else if (vh.scaf[i].type(vh, i) === 'end') {
                if (vh.stap[i].type(vh, i) === 'empty') {
                    nodes.add_end(i);
                } else if (vh.stap[i].type(vh, i) === 'continue') {
                    nodes.add_begin(i + 1 * direction);
                    nodes.add_end(i);
                } else if (vh.stap[i].type(vh, i) === 'begin') {
                    nodes.add_end(i);
                } else if (vh.stap[i].type(vh, i) === 'end') {
                    nodes.add_begin(i + 1 * direction);
                    nodes.add_end(i - 1 * direction);
                }
            } else if (vh.scaf[i].type(vh, i) === 'continue') {
                if (vh.stap[i].type(vh, i) === 'begin') {
                    nodes.add_begin(i + 1 * direction);
                    nodes.add_end(i);
                } else if (vh.stap[i].type(vh, i) === 'end') {
                    nodes.add_begin(i);
                    nodes.add_end(i - 1 * direction);
                }
            }
        }
    }
    return nodes;
}


function generate_vhelices_origami_sq(vhelix_direction, vhelix_perp, h) {
    let g = new cu.StrandGenerator();
    //  generate helix angles
    let helix_angles: number[] = new Array(h.len - 1);
    //  hard upper limit on pitch angle seems to be between 54.5 and 55 degrees
    for(let i=0; i<helix_angles.length; i++) {
        let modi = i % 32;
        if (modi < 2) {
            helix_angles[i] = 28 * Math.PI / 180;
        } else if (modi === 2) {
            helix_angles[i] = 36 * Math.PI / 180;
        } else if (modi === 3) {
            helix_angles[i] = 54.375 * Math.PI / 180;
        } else if (modi === 4) {
            helix_angles[i] = 37 * Math.PI / 180;
        } else if ([5, 6].includes(modi)) {
            helix_angles[i] = 27.6666666666666 * Math.PI / 180;
        } else if (modi === 7) {
            helix_angles[i] = 30.6666666666666 * Math.PI / 180;
        } else if ([8, 9].includes(modi)) {
            helix_angles[i] = 29.3333333333 * Math.PI / 180;
        } else if (modi === 10) {
            helix_angles[i] = 34.3333333333 * Math.PI / 180;
        } else if (modi === 11) {
            helix_angles[i] = 54.5 * Math.PI / 180;
        } else if ([12, 13].includes(modi)) {
            helix_angles[i] = (28.91666666666 * Math.PI / 180);
        } else if ([14, 15, 16, 17].includes(modi)) {
            helix_angles[i] = 31.16666666666 * Math.PI / 180;
        } else if (modi === 18) {
            helix_angles[i] = 35.5 * Math.PI / 180;
        } else if (modi === 19) {
            helix_angles[i] = 52 * Math.PI / 180;
        } else if (modi === 20) {
            helix_angles[i] = 35.5 * Math.PI / 180;
        } else if ([21, 22].includes(modi)) {
            helix_angles[i] = 27.5 * Math.PI / 180;
        } else if (modi === 23) {
            helix_angles[i] = 35.5 * Math.PI / 180;
        } else if (modi >= 24 && modi < 27) {
            helix_angles[i] = 30 * Math.PI / 180;
        } else if (modi === 27) {
            helix_angles[i] = 52 * Math.PI / 180;
        } else if (modi === 28) {
            helix_angles[i] = 35.5 * Math.PI / 180;
        } else {
            helix_angles[i] = 30.91666666666 * (Math.PI / 180);
        }
    }

    //  make sure the helices are periodic in 32 bases
    let total_sum = 0;
    for(let i=0; i<31; i++) {
        total_sum += helix_angles[i];
    }

    for(let i=0; i<helix_angles.length; i++) {
        if (i % 32 === 31) {
            helix_angles[i] = 1080 * Math.PI / 180 - total_sum;
        }
    }

    //  make the virtual helices
    let pos: THREE.Vector3, direction: THREE.Vector3, perp: THREE.Vector3, rot: number, angles: number[];
    if (h.num % 2 === 0) {
        pos = new THREE.Vector3(h.col * DIST_SQUARE, h.row * DIST_SQUARE, 0);
        direction = vhelix_direction.clone();
        perp = vhelix_perp.clone();
        rot = 0.;
        angles = helix_angles.slice();
    } else {
        pos = new THREE.Vector3(h.col * DIST_SQUARE, h.row * DIST_SQUARE, (h.len - 1) * base.BASE_BASE);
        direction = vhelix_direction.clone().negate();
        perp = vhelix_perp.clone().negate();
        rot = -(utils.sum(...helix_angles)) % (2 * Math.PI);
        angles = helix_angles.slice().reverse();
    }

    const strands = g.generate_or_sq(h.len, pos, direction, perp, true, rot, angles)

    return [strands, helix_angles, pos, rot, direction, perp];
}


function generate_vhelices_origami_he(vhelix_direction: THREE.Vector3, vhelix_perp: THREE.Vector3, h: vhelix) {
    let g = new cu.StrandGenerator();
    //  generate helix angles
    let helix_angles = new Array(h.len - 1);

    for(let i=0; i<helix_angles.length; i++) {
        const modi = i % 21;
        if (modi === 0) {
            helix_angles[i] = 32.571 * Math.PI / 180;
        } else if (modi === 1) {
            helix_angles[i] = 36 * Math.PI / 180;
        } else if ([1, 2, 3].includes(modi)) {
            helix_angles[i] = 42 * Math.PI / 180;
        } else if ([5, 6, 7].includes(modi)) {
            helix_angles[i] = 29.143 * Math.PI / 180;
        } else if (modi === 8) {
            helix_angles[i] = 32 * Math.PI / 180;
        } else if ([9, 10].includes(modi)) {
            helix_angles[i] = 44 * Math.PI / 180;
        } else if ([12, 13, 14].includes(modi)) {
            helix_angles[i] = 28.571 * Math.PI / 180;
        } else if ([16, 17].includes(modi)) {
            helix_angles[i] = 41.5 * Math.PI / 180;
        } else if ([19, 20].includes(modi)) {
            helix_angles[i] = 28.476 * Math.PI / 180;
        } else {
            helix_angles[i] = 720. / 21 * (Math.PI / 180.);
        }
    }

    //  make sure it's periodic
    let total_sum = 0;
    for(let i=0; i<20; i++) {
        total_sum += helix_angles[i];
    }

    for(let i=0; i<helix_angles.length; i++) {
        if (i % 21 === 20) {
            helix_angles[i] = 720. * Math.PI / 180 - total_sum;
        }
    }

    //  make the virtual helices
    let pos: THREE.Vector3, direction: THREE.Vector3, perp: THREE.Vector3, rot: number, strands;
    if (h.num % 2 === 0) {
        pos = new THREE.Vector3(h.col * Math.sqrt(3) * DIST_HEXAGONAL / 2, h.row * 3 * DIST_HEXAGONAL / 2, 0);
        direction = vhelix_direction.clone();
        perp = vhelix_perp.clone();
        rot = 0.;
        strands = g.generate_or_sq(h.len, pos, direction, perp, true, rot, helix_angles);
    } else {
        pos = new THREE.Vector3(h.col * Math.sqrt(3) * DIST_HEXAGONAL / 2, h.row * 3 * DIST_HEXAGONAL / 2 + DIST_HEXAGONAL / 2, (h.len - 1) * base.BASE_BASE);
        direction = vhelix_direction.clone().negate();
        perp = vhelix_perp.clone().negate();
        if (base.MM_GROOVING) {
            rot = -utils.sum(...helix_angles) % (2 * Math.PI) - 0.07;
        } else {
            rot = -utils.sum(...helix_angles) % (2 * Math.PI);
        }
        let angles = helix_angles.slice().reverse();
        strands = g.generate_or_sq(h.len, pos, direction, perp, true, rot, angles);
    }

    return [strands, helix_angles, pos, rot, direction, perp];
}


//  cadnano object structure
class vstrands {
    vhelices: vhelix[];

    constructor() {
        this.vhelices = [];
    }

    add_vhelix(toadd) {
        this.vhelices.push(toadd);
    }

    bbox() {
        let rows = [];
        let cols = [];
        let lens = [];

        for (const h of this.vhelices) {
            rows.push(h.row);
            cols.push(h.col);
            lens.push(h.stap.length);
        }

        const dr = DIST_SQUARE * (Math.max(...rows) - Math.min(...rows) + 2);
        const dc = DIST_SQUARE * (Math.max(...cols) - Math.min(...cols) + 2);
        const dl = 0.34 * (Math.max(...lens) + 2);
        
        return 2 * Math.max(dr, dc, dl) * BOX_FACTOR;
    }
    
    toString() {
        let a = '{\n"vstrands":[\n';
        if (this.vhelices.length > 0) {
            for (const h of this.vhelices) {
                a += `${h},`;
            }
            a = a.slice(0, a.length - 1);
        }
        a += '}\n';
        return a;
    }
}


class vhelix {
    stapLoop: number[];
    scafLoop: number[];
    skip: number[];
    loop: number[];
    stap_colors: [number, number][];
    row: number;
    col: number;
    num: number;
    stap: square[];
    scaf: square[];
    cad_index: number;
    skiploop_bases: number;
    constructor() {
        this.stapLoop = [];
        this.scafLoop = [];
        this.skip = [];
        this.loop = [];
        this.stap_colors = [];
        this.row = 0;
        this.col = 0;
        this.num = 0;
        this.stap = [];
        this.scaf = [];
        this.cad_index = -1;
        this.skiploop_bases = 0;
    }

    get len() {
        return Math.max(this.scaf.length, this.stap.length);
    }

    add_square(toadd, which) {
        if (which === 'stap') {
            this.stap.push(toadd);
        } else if (which === 'scaf') {
            this.scaf.push (toadd);
        } else {
            base.Logger.log("Cannot add square that is not scaf or stap. Dying now", base.Logger.CRITICAL);
            //sys.exit(1);
        }
    }

    toString() {
        let a = '{\n';
        a += '"stapLoop":[';
        if (this.stapLoop.length > 0) {
            for (const i of this.stapLoop) {
                a += `${i},`;
            }
            a = a.slice(0, a.length - 1)  //  remove last comma;
        }
        a += '],\n';
        a += '"skip":[';
        if (this.skip.length > 0) {
            for (const e of this.skip) {
                a += `${e},`;
            }
            a = a.slice(0, a.length - 1)  //  remove last comma;
        }
        a += '],\n';
        a += '"loop":[';
    
        if (this.loop.length > 0) {
            for (const e of this.loop) {
                a += `${e},`;
            }
            a = a.slice(0, a.length - 1)  //  remove last comma;
        }
        a += '],\n';
        a += '"stap_colors":[';
        if (this.stap_colors.length > 0) {
            for (const e of this.stap_colors) {
                a += `${e},`;
            }
            a = a.slice(0, a.length - 1);  //  remove last comma;
        }
        a += '],\n';

        a += `"row":${this.row},\n`;
        a += `"col":${this.col},\n`;
        a += `"num":${this.num},\n`;
        
        a += '"scafLoop":[';
        if (this.scafLoop.length > 0) {
            for (const i of this.scafLoop) {
                a += `${i},`;
            }
            a = a.slice(0, a.length - 1)  //  remove last comma;
        }
        a += '],\n';
        
        a += '"stap":[';
        if (this.stap.length > 0) {
            for (const i of this.stap) {
                a += `${i},`;
            }
            a = a.slice(0, a.length - 1)  //  remove last comma;
        }
        a += '],\n';
        
        a += '"scaf":[';
        if (this.scaf.length > 0) {
            for (const i of this.scaf) {
                a += `${i},`;
            }
            a = a.slice(0, a.length - 1)  //  remove last comma;
        }
        a += ']\n}';
        return a;
    }
}


class square {
    V_0: number;
    b_0: number;
    V_1: number;
    b_1: number;
    constructor(V_0=-1, b_0=-1, V_1=-1, b_1=-1) {
        /*
        V_0, b_0, V_1, b_1 are integer indices correspond to:
        virtual_helix_behind, virtual_base_behind, virtual_helix_ahead, virtual_base_ahead
        */
        this.V_0 = V_0;
        this.b_0 = b_0;
        this.V_1 = V_1;
        this.b_1 = b_1;
    }

    toString() {
        return `[${this.V_0},${this.b_0},${this.V_1},${this.b_1}]`;
    }

    type(vhelix: vhelix, myid) {
        //  find type of strand (junction) on this square
        //  currently direction always equals zero...
        let direction = 0;
        if (this.V_0 === -1 && this.b_0 === -1) {
            if (this.V_1 === -1 && this.b_1 === -1) {
                return 'empty';
            } else if (this.V_1 === vhelix.num && Math.abs(this.b_1 - myid) === 1) {
                if (direction === 0) {
                    return 'begin';
                } else {
                    return 'end';
                }
            }
        } else if (this.V_0 === vhelix.num && Math.abs(this.b_0 - myid) === 1) {
            if (this.V_1 === -1) {
                if (direction === 0) {
                    return 'end';
                } else {
                    return 'begin';
                }
            } else if (this.V_1 === vhelix.num && Math.abs(this.b_1 - myid) === 1) {
                return 'continue'
            } else {
                //  join
                if (direction === 0) {
                    return 'end';
                } else {
                    return 'begin';
                }
            }
        } else {
            if (this.V_1 === vhelix.num && Math.abs(this.b_1 - myid) === 1) {
                if (direction === 0) {
                    return 'begin'
                } else {
                    return 'end'
                }
            }
        }
        //  shouldn't get to here
        base.Logger.log('unexpected square array', base.Logger.WARNING);
    }
}

function parse_cadnano(json_string) {
    let cadsys = new vstrands();
    
    let cadnano = JSON.parse(json_string);
    for (const vstrand of cadnano["vstrands"]) {
        let vh = new vhelix();
        for (const [key, val] of Object.entries(vstrand)) {
            if (key === "skip") {
                vh.skip = (val as number[]).map(x=>Math.abs(x));
            } else if (key === "stap" || key === "scaf") {
                vh[key] = (val as number[][]).map(x=>new square(...x));
            } else {
                vh[key] = val;
            }
        }
        vh.skiploop_bases = vh.skip.length + utils.sum(...vh.loop) - utils.sum(...vh.skip);
        cadsys.add_vhelix(vh);
    }
    return cadsys;
}

// https://stackoverflow.com/a/35891015
function getMostCommon(array) {
    var count = {};
    array.forEach(function (a) {
        count[a] = (count[a] || 0) + 1;
    });
    return Object.keys(count).reduce(function (r, k, i) {
        if (!i || count[k] > count[r[0]]) {
            return [k];
        }
        if (count[k] === count[r[0]]) {
            r.push(k);
        }
        return r;
    }, []);
}

function loadCadnano(source_file: string, grid: string, scaffold_seq?: string, side?: number) {
    let origami_sq = false;
    let origami_he = false;
    if (grid === "sq") {
        origami_sq = true;
    } else if (grid === "he") {
        origami_he = true;
    } else {
        alert("Lattice_type should be either 'sq' or 'he'");
        return;
    }

    let vh_vb2nuc = new cu.vhelix_vbase_to_nucleotide();
    let vh_vb2nuc_final = new cu.vhelix_vbase_to_nucleotide();

    let cadsys = parse_cadnano(source_file);

    let vhelix_counter = 0;
    if (side === undefined) {
        side = cadsys.bbox();
        base.Logger.log(`Using default box size, a factor ${BOX_FACTOR} larger than the size of the cadnano system`, base.Logger.INFO);
    }
    let vhelix_direction_initial = new THREE.Vector3(0, 0, 1);
    let vhelix_perp_initial = new THREE.Vector3(1, 0, 0);
    if (origami_sq) {
        vhelix_perp_initial = vhelix_rotation_origami_sq(vhelix_direction_initial, vhelix_perp_initial);
    } else if (origami_he) {
        vhelix_perp_initial = vhelix_rotation_origami_he(vhelix_direction_initial, vhelix_perp_initial);
    }

    let slice_sys = new base.System(new THREE.Vector3(side, side, side));
    let final_sys = new base.System(new THREE.Vector3(side, side, side));
    let strand_number = -1;
    let partner_list_scaf = [];
    let partner_list_stap = [];
    let found_partner = false;
    let join_list_scaf = [];
    let join_list_stap = [];
    let begin_helix = -1;
    let end_helix = -1;
    for (const h of cadsys.vhelices) {
        h.cad_index = vhelix_counter;
        let strands, helix_angles, pos, rot, vhelix_direction, vhelix_perp;
        if (origami_sq) {
            [strands, helix_angles, pos, rot, vhelix_direction, vhelix_perp] = generate_vhelices_origami_sq(vhelix_direction_initial, vhelix_perp_initial, h);
        } else if (origami_he) {
            [strands, helix_angles, pos, rot, vhelix_direction, vhelix_perp] = generate_vhelices_origami_he(vhelix_direction_initial, vhelix_perp_initial, h);
        }

        let nodes = build_nodes(h);

        //  read the scaffold squares and add strands to slice_sys
        let i = 0;
        for (const s of h.scaf) {
            if (s.V_0 === -1 && s.b_0 === -1) {
                if (s.V_1 === -1 && s.b_0 === -1) {
                    ;
                } else if (s.V_1 === h.num && Math.abs(s.b_1 - i) === 1) {
                    if (h.num % 2 === 0) {
                        strand_number++;
                    }
                    begin_helix = i;
                    if (h.num % 2 === 1) {
                        slice_sys = add_slice(slice_sys, h, begin_helix, end_helix, nodes, strands, pos, vhelix_direction, vhelix_perp, rot, helix_angles, 0);
                        vh_vb2nuc = add_slice_nupack(h, strand_number, begin_helix, end_helix, vh_vb2nuc, 2);
                    }
                } else {
                    base.Logger.log("unexpected square array", base.Logger.WARNING);
                }
            } else if (s.V_0 === h.num && Math.abs(s.b_0 - i) === 1) {
                if (s.V_1 === -1 && s.b_1 === -1) {
                    if (h.num % 2 === 1) {
                        strand_number++;
                    }
                    end_helix = i;
                    if (h.num % 2 === 0) {
                        slice_sys = add_slice(slice_sys, h, begin_helix, end_helix, nodes, strands, pos, vhelix_direction, vhelix_perp, rot, helix_angles, 0);
                        vh_vb2nuc = add_slice_nupack(h, strand_number, begin_helix, end_helix, vh_vb2nuc, 2);
                    }
                } else if (s.V_1 === h.num && Math.abs(s.b_1 - i) === 1) {
                    ;
                } else {
                    if (h.num % 2 === 1) {
                        strand_number++;
                    }
                    end_helix = i;
                    if (h.num % 2 === 0 ) {
                        slice_sys = add_slice(slice_sys, h, begin_helix, end_helix, nodes, strands, pos, vhelix_direction, vhelix_perp, rot, helix_angles, 0);
                        vh_vb2nuc = add_slice_nupack(h, strand_number, begin_helix, end_helix, vh_vb2nuc, 2);
                    }
                    let column;
                    if (h.num % 2 === 1) {
                        column = i;
                    } else {
                        column = i;
                    }
                    for (let j=0; j<partner_list_scaf.length; j++) {
                        if (utils.arraysEqual([h.num, column], partner_list_scaf[j])) {
                            join_list_scaf[j].splice(0, 0, strand_number);
                            found_partner = true;
                        }
                    }
                    if (!found_partner ) {
                        join_list_scaf.push([strand_number])
                        partner_list_scaf.push([s.V_1, s.b_1])
                    }
                    found_partner = false;
                }
            } else {
                if (s.V_1 === -1 && s.b_1 === -1) {
                    base.Logger.log("unexpected square array", base.Logger.WARNING)
                } else if (s.V_1 === h.num && Math.abs(s.b_1 - i) === 1) {
                    if (h.num % 2 === 0) {
                        strand_number++;
                    }
                    begin_helix = i;
                    if (h.num % 2 === 1) {
                        slice_sys = add_slice(slice_sys, h, begin_helix, end_helix, nodes, strands, pos, vhelix_direction, vhelix_perp, rot, helix_angles, 0);
                        vh_vb2nuc = add_slice_nupack(h, strand_number, begin_helix, end_helix, vh_vb2nuc, 2);
                    }
                    for (let j=0; j<partner_list_scaf.length; j++) {
                        let column;
                        if (h.num % 2 === 1) {
                            column = i;
                        } else {
                            column = i;
                        }
                        if (utils.arraysEqual([h.num, column], partner_list_scaf[j])) {
                            join_list_scaf[j].push(strand_number)
                            found_partner = true;
                        }
                    }
                    if (!found_partner ) {
                        join_list_scaf.push([strand_number])
                        partner_list_scaf.push([s.V_0, s.b_0])
                    }
                    found_partner = false;
                } else {
                    base.Logger.log("unexpected square array", base.Logger.WARNING);
                }
            }
            i++;
        }
        if (slice_sys.N_strands === 0) {
            base.Logger.log(`No scaffold strand found in virtual helix n. ${h.num}: staples-only virtual helices are not supported`, base.Logger.WARNING)
            continue;
        }

        //  read the staple squares and add strands to slice_sys
        i = 0;
        for (const s of h.stap) {
            if (s.V_0 === -1 && s.b_0 === -1) {
                if (s.V_1 === -1 && s.b_0 === -1) {
                    ;
                } else if (s.V_1 === h.num && Math.abs(s.b_1 - i) === 1) {
                    if (h.num % 2 === 1) {
                        strand_number++;
                    }
                    begin_helix = i;
                    if (h.num % 2 === 0) {
                        slice_sys = add_slice(slice_sys, h, begin_helix, end_helix, nodes, strands, pos, vhelix_direction, vhelix_perp, rot, helix_angles, 1);
                        vh_vb2nuc = add_slice_nupack(h, strand_number, begin_helix, end_helix, vh_vb2nuc, 3);
                    }
                } else {
                    base.Logger.log("unexpected square array", base.Logger.WARNING);
                }
            } else if (s.V_0 === h.num && Math.abs(s.b_0 - i) === 1) {
                if (s.V_1 === -1 && s.b_1 === -1) {
                    if (h.num % 2 === 0) {
                        strand_number++;
                    }
                    end_helix = i;
                    if (h.num % 2 === 1) {
                        slice_sys = add_slice(slice_sys, h, begin_helix, end_helix, nodes, strands, pos, vhelix_direction, vhelix_perp, rot, helix_angles, 1);
                        vh_vb2nuc = add_slice_nupack(h, strand_number, begin_helix, end_helix, vh_vb2nuc, 3);
                    }
                } else if (s.V_1 === h.num && Math.abs(s.b_1 - i) === 1) {
                    ;
                } else {
                    if (h.num % 2 === 0) {
                        strand_number++;
                    }
                    end_helix = i;
                    if (h.num % 2 === 1) {
                        slice_sys = add_slice(slice_sys, h, begin_helix, end_helix, nodes, strands, pos, vhelix_direction, vhelix_perp, rot, helix_angles, 1);
                        vh_vb2nuc = add_slice_nupack(h, strand_number, begin_helix, end_helix, vh_vb2nuc, 3);
                    }
                    let column = i;
                    for (let j=0; j<partner_list_stap.length; j++) {
                        if (utils.arraysEqual([h.num, column], partner_list_stap[j])) {
                            join_list_stap[j].splice(0, 0, strand_number);
                            found_partner = true;
                        }
                    }
                    if (!found_partner ) {
                        join_list_stap.push([strand_number]);
                        partner_list_stap.push([s.V_1, s.b_1]);
                    }
                    found_partner = false;
                }
            } else {
                if (s.V_1 === -1 && s.b_1 === -1) {
                    base.Logger.log("unexpected square array", base.Logger.WARNING)
                } else if (s.V_1 === h.num && Math.abs(s.b_1 - i) === 1) {
                    if (h.num % 2 === 1) {
                        strand_number++;
                    }
                    begin_helix = i;
                    if (h.num % 2 === 0) {
                        slice_sys = add_slice(slice_sys, h, begin_helix, end_helix, nodes, strands, pos, vhelix_direction, vhelix_perp, rot, helix_angles, 1);
                        vh_vb2nuc = add_slice_nupack(h, strand_number, begin_helix, end_helix, vh_vb2nuc, 3);
                    }
                    for (let j=0; j<partner_list_stap.length; j++) {
                        let column = i;
                        if (utils.arraysEqual([h.num, column], partner_list_stap[j])) {
                            join_list_stap[j].push(strand_number)
                            found_partner = true;
                        }
                    }
                    if (!found_partner ) {
                        join_list_stap.push([strand_number])
                        partner_list_stap.push([s.V_0, s.b_0])
                    }
                    found_partner = false;
                } else {
                    base.Logger.log("unexpected square array", base.Logger.WARNING);
                }
            }
            i++;
        }
        vhelix_counter++;
    }

    let join_lists = [join_list_scaf, join_list_stap];

    //  add strands to final_sys that aren't joined
    let join_list_unpacked = [];
    for (let a=0; a<2; a++) {
        for (const i of join_lists[a]) {
            join_list_unpacked.push(...i);
        }
    }
    for(let i=0; i<slice_sys._strands.length; i++) {
        if (!join_list_unpacked.includes(i)) {
            final_sys.add_strand(slice_sys._strands[i]);
            vh_vb2nuc_final.add_strand(i, vh_vb2nuc);
        }
    }

    let restart: boolean;

    for (let a=0; a<2; a++) {
        const join_list = join_lists[a];
        let all_are_joined = false;
        restart = false;

        //  check distance between the backbones we are about to join
        for (const pair of join_list) {
            let strand1 = slice_sys._strands[pair[0]];
            let strand2 = slice_sys._strands[pair[1]];
            let backbone_backbone_dist = strand1._nucleotides.slice(-1)[0].distance(strand2._nucleotides[0], false);
            let absolute_bb_dist = Math.sqrt(backbone_backbone_dist.dot(backbone_backbone_dist));
            if (absolute_bb_dist > 1.0018 || absolute_bb_dist < 0.5525) {
                base.Logger.log(`The backbone-backbone distance across joints is ${absolute_bb_dist}: it will have to be relaxed with preliminary simulations`, base.Logger.WARNING);
            }
        }

        //  match up all the pairs of joins that involve the same strand
        let circular = [];
        while (!all_are_joined ) {
            restart = false;
            for(let i=0; i<join_list.length; i++) {
                if (restart) {
                    break;
                }
                for (let j=0; j<join_list.length; j++) {
                    if (restart) {
                        break;
                    }
                    if (join_list[i][0] === join_list[j].slice(-1)[0]) {
                        if (i != j) {
                            join_list[j].push(...join_list[i].slice(1, ))
                            join_list.splice(i,1);
                            restart = true;
                            break;
                        } else {
                            if (!circular.includes(i)) {
                                circular.push(i)
                            }
                        }
                    }
                }
            }

            if (!restart ) {
                all_are_joined = true;
            }
        }

        //  add joined strands
        for (let ii=0; ii<join_list.length; ii++) {
            const join = join_list[ii];
            let joined_strand = slice_sys._strands[join[0]];
            if (circular.includes(ii)) {
                for (const k of utils.range(1, join.length - 1)) {
                    joined_strand = joined_strand.append(slice_sys._strands[join[k]]);
                }
                joined_strand.make_circular(true);
            } else {
                for (const k of utils.range(1, join.length)) {
                    joined_strand = joined_strand.append(slice_sys._strands[join[k]]);
                }
            }

            final_sys.add_strand(joined_strand);

            //  This is a bug fix. Ben 12/2/14
            //  for a circular strand we need to terminate the strand one element early (so reduce the length
            //  of the range by 1), since the final element is just a repeat of the first one.
            let joining_range: number[];
            if (joined_strand._circular) {
                joining_range = utils.range(join.length - 2);
            } else {
                joining_range = utils.range(join.length - 1);
            }
            //  add joined strands to v2n index
            for (const k of joining_range) {
                vh_vb2nuc_final.add_strand(join[k], vh_vb2nuc, true);
            }
            vh_vb2nuc_final.add_strand(join[joining_range[joining_range.length-1] + 1], vh_vb2nuc, false);
        }
    }

    //  Fix to reverse the direction of every strand so that the 3' to 5' direction is the same
    //  as in Cadnano. In cadnano the strands point in the 5' to 3' direction, whereas in oxDNA
    //  they point in the 3' to 5' direction. Ben 29/11/13
    let rev_sys = new base.System(final_sys._box);
    for (const strand of final_sys._strands) {
        let reverse_nucs = strand._nucleotides.map(nuc=>nuc);
        reverse_nucs.reverse()
        let rev_strand = new base.Strand();
        for (const nuc of reverse_nucs) {
            rev_strand.add_nucleotide(new base.Nucleotide(nuc.cm_pos, nuc._a1, nuc._a3.clone().negate(), nuc._base));
        }
        if (strand._circular) {
            rev_strand.make_circular();
        }
        rev_sys.add_strand(rev_strand);
    }
    //  also reverse the vhelix_vbase_to_nucleotide order so it corresponds to the reversed system
    let vh_vb2nuc_rev = new cu.vhelix_vbase_to_nucleotide();
    //  count the number of nucleotides up to but not including the nucleotides in strand ii
    let nnucs_to_here = utils.range(rev_sys._N_strands);
    let nuc_total = 0;
    for (let strandii=0; strandii<rev_sys._strands.length; strandii++) {
        const strand = rev_sys._strands[strandii];
        nnucs_to_here[strandii] = nuc_total;
        nuc_total += strand._nucleotides.length;
    }

    let id_to_pos = new Map();
    //  fill in the _scaf and _stap dicts for the reverse vhelix_vbase_to_nucleotide object

    for (const [vh, vb] of vh_vb2nuc_final._scaf.keys()) {
        let [strandii, nuciis] = vh_vb2nuc_final._scaf.get([vh, vb]);
        let rev_nuciis = [];
        for (const nucii of nuciis) {
            const nuc = rev_sys._strands[strandii]._nucleotides.length - 1 - (nucii - nnucs_to_here[strandii]) + nnucs_to_here[strandii];
            rev_nuciis.push(nuc)
            id_to_pos.set(nuc, [vh, vb]);
        }
        vh_vb2nuc_rev.add_scaf(vh, vb, strandii, rev_nuciis);
    }
    for (const [vh, vb] of vh_vb2nuc_final._stap.keys()) {
        let [strandii, nuciis] = vh_vb2nuc_final._stap.get([vh, vb]);
        let rev_nuciis = [];
        for (const nucii of nuciis) {
            const nuc = rev_sys._strands[strandii]._nucleotides.length - 1 - (nucii - nnucs_to_here[strandii]) + nnucs_to_here[strandii];
            rev_nuciis.push(nuc)
            id_to_pos.set(nuc, [vh, vb]);
        }
        vh_vb2nuc_rev.add_stap(vh, vb, strandii, rev_nuciis);
    }

    if (rev_sys.N === 0) {
        base.Logger.log("The generated configuration is empty: this might be due to this conversion module not supporting virtual helices containing no scaffold strands.", base.Logger.CRITICAL)
        return;
    }

    //  Find colors
    let pos_to_color = new cu.PairMap();
    for (const vh of cadsys.vhelices) {
        for (const [vb, color] of vh.stap_colors) {
            pos_to_color.set([vh.num, vb], color);
        }
    }

    //  Create inverse mapping to find pairs by their positions
    let pos_to_id = new cu.PairMap();
    for (const [bid, helix] of id_to_pos.entries()) {
        if (!pos_to_id.has(helix)) {
            pos_to_id.set(helix, []);
        }
        pos_to_id.get(helix).push(bid);
    }

    //  Create mapping to find nucleotide by its index
    let id_to_nucleotide = new Map()
    for(const s of rev_sys._strands) {
        for(const n of s._nucleotides) {
            id_to_nucleotide.set(n.index, n);
        }
    }

    //  Find the index of the scaffold strand (there really should be some
    //  easier way of doing this)
    let strands = [];
    for (const nucIds of pos_to_id.values()) {
        for(const nucId of nucIds) {
            try {
                strands.push(id_to_nucleotide.get(nucId).strand);
            } catch (error) {
                console.log(`Could not find nucId ${nucId} strand`);
            }
        }
    }
    let scaffold_index = getMostCommon(strands)[0];

    //  Go through system and add extra information
    //  (basepairs, clusters and colors)
    for (const s of rev_sys._strands) {
        for (const n of s._nucleotides) {
            if (id_to_pos.has(n.index)) {
                //  Find which helix and position the nucleotide had
                let [vh, vb] = id_to_pos.get(n.index);
                let paired = pos_to_id.get([vh,vb]);

                //  If double-stranded, save basepair
                if (paired.length > 1) {
                    n.pair = id_to_nucleotide.get(paired.filter(idx=>idx != n.index)[0]);
                }

                //  Get the staple strand nucleotide at this position
                let staple_ids = pos_to_id.get([vh,vb]).filter(n=>(id_to_nucleotide.get(n).strand != scaffold_index));
                if (staple_ids.length > 0) {
                    let staple_nuc: base.Nucleotide = id_to_nucleotide.get(staple_ids[0]);

                    //  If this is the staple strand, check if it has a color
                    if (n === staple_nuc) {
                        if (pos_to_color.has([vh, vb])) {
                            n.color = pos_to_color.get([vh, vb]);
                        }
                    } else {
                        n.color = 3633362; // Make scaffold blue
                    }
                }
            }
        }
    }
    //  If there's any colored nucleotide in a strand,
    //  use that color for the whole strand
    for (const s of rev_sys._strands) {
        for (const n of s._nucleotides) {
            if (n.color !== undefined) {
                for (const other of s._nucleotides) {
                    other.color = n.color;
                }
                break;
            }
        }
    }

    rev_sys.calc_clusters();

    // Set scaffold sequence
    let scaffold_strand = rev_sys._strands[scaffold_index];
    let scaffold_bases: number[];
    if (scaffold_seq) {
        if (scaffold_strand.N > scaffold_seq.length) {
            base.Logger.log(
                `Provided scaffold sequence is ${scaffold_seq.length}nt `+
                `but needs to be at least ${scaffold_strand.N} to cover the scaffold. `
            );
        } else {
            base.Logger.log("Applying custom sequence");
            scaffold_bases = Array.from(scaffold_seq).map(b=>base.base_to_number[b]);
        }
    }
    if (scaffold_bases === undefined) {
        base.Logger.log("Applying random sequence");
        scaffold_bases = utils.randint(0, 4, scaffold_strand.N) as number[];
    }

    for (let nuc_i=0; nuc_i<scaffold_strand.N; nuc_i++) {
        let n = scaffold_strand._nucleotides[nuc_i];
        n._base = scaffold_bases[scaffold_strand.N-nuc_i-1];
        if (n.pair) {
            n.pair._base = 3 - n._base;
        }
    }

    return rev_sys
}

export {loadCadnano}