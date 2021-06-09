import * as THREE from "three";
declare const FROM_ANGSTROM_TO_OXDNA: number;
declare class Nucleotide {
    base: string;
    idx: number;
    name: string;
    base_atoms: Atom[];
    phosphate_atoms: Atom[];
    sugar_atoms: Atom[];
    named_atoms: {};
    ring_names: string[];
    chain_id: string;
    a1: THREE.Vector3;
    a3: THREE.Vector3;
    a2: THREE.Vector3;
    check: number;
    constructor(name: string, idx: number);
    get atoms(): any[];
    add_atom(a: any): void;
    get_com(atoms?: any): THREE.Vector3;
    compute_a3(): void;
    compute_a1(): void;
    compute_as(): void;
    correct_for_large_boxes(box: any): void;
    to_pdb(chain_identifier: any, print_H: any, residue_serial: any, residue_suffix: any, residue_type: any): string;
    to_mgl(): string;
    rotate(R: any): void;
    set_com(new_com: any): void;
}
declare class Atom {
    name: string;
    alternate: string;
    residue: string;
    chain_id: string;
    residue_idx: number;
    pos: THREE.Vector3;
    static serial_atom: number;
    constructor(pdb_line: string);
    clone(): Atom;
    is_hydrogen(): void;
    shift(diff: THREE.Vector3): void;
    to_pdb(chain_identifier: any, residue_serial: any, residue_suffix: any): string;
}
export { Nucleotide, Atom, FROM_ANGSTROM_TO_OXDNA };
