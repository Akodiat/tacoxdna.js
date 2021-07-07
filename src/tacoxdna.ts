import {loadCadnano} from './cadnanoLoader'
import {loadRpoly} from './rpolyLoader'
import {loadTiamat} from './tiamatLoader'
import {loadPDB} from './pdbLoader'
import {Logger, System} from './libs/base'

function convertFromTo(inputs: string[], from: string, to: string, opts) {
    let sys: System;
    switch (from) {
        case 'cadnano':
            sys = loadCadnano(inputs[0], opts.grid, opts.sequence, opts.side);
            break;
        case 'rpoly':
            sys = loadRpoly(inputs[0], opts.sequence);
            break;
        case 'tiamat':
            sys = loadTiamat(inputs[0], opts.tiamat_version, opts.isDNA, opts.default_val);
            break;
        case 'pdb':
            sys = loadPDB(inputs[0], opts.strand_dir, opts.models_as_strands);
            break;
        default:
            console.log("Unknown input format: "+from);
            return;
    }
    switch (to) {
        case 'oxview':
            return sys.print_oxview_output();
        case 'oxdna':
            return sys.print_lorenzo_output();
        default:
            console.log("Unknown output format: "+from);
            break;
    }
}

export {convertFromTo, Logger}