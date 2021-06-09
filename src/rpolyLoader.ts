import * as THREE from 'three'
import * as base from "./libs/base";
import * as cu from "./libs/cadnano_utils";


// rpoly file contains center coordinate of helix,
// "generate" needs end coordiates of helix
function move_along_vector(point: THREE.Vector3, vector: THREE.Vector3, length) {
    // 0.4 is the length of a base pair in oxDNA units, move half the helixlength down
    let move_distance = parseFloat(length) * 0.4 / 2.0;
    return point.clone().sub(vector.clone().multiplyScalar(move_distance));
}

function loadRpoly(source_file: string) {
    // Read File
    // 'data' stores helix coordinates + rotaion in quaternion
    let data = [];
    
    let rev_helix_connections = [];  // staple connection information,
    let fwd_helix_connections = [];  // scaffold connections
    let count = 0;
    
    for(const line of source_file.split('\n')) {
        if (line.startsWith('hb')) {
            data.splice(count, 0, line.split(' '));
            count += 1;
        } else if (line.startsWith('c')) {
            if (!line.includes('f3')) {
                rev_helix_connections.push([
                    parseInt(line.match(/c helix_(.+?) /)[1]),
                    parseInt(line.match(/\' helix_(.+?) /)[1])
                ]);  // Extract connection information
            } else {
                fwd_helix_connections.push([
                    parseInt(line.match(/c helix_(.+?) /)[1]),
                    parseInt(line.match(/\' helix_(.+?) /)[1])
                ]);
            }
        }
    }
    
    let generator = new cu.StrandGenerator();
    
    // temporary system to store staple fragments before later connecting them
    let staple_fragments = new base.System(new THREE.Vector3(100, 100, 100));
    let scaffold_fragments = new base.System(new THREE.Vector3(100,100,100))
    
    // Reads orientation from the "data" and produces rotations from the Quaternian coordinates
    let largest_size = 0;
    for (let n=0; n<data.length; n++) {
        const i = data[n];
    
        // 0.84 scaling is ad hoc solution to get good looking models
        let position = new THREE.Vector3(
            parseFloat(i[3]) / 0.84,
            parseFloat(i[4]) / 0.84,
            parseFloat(i[5]) / 0.84
        );
        
        let n_bp = parseInt(i[2])
        
        let q = new THREE.Quaternion(parseFloat(i[6]), parseFloat(i[7]), parseFloat(i[8]), parseFloat(i[9]))  // find the helix roation Info from file
        let vec = new THREE.Vector3(0.0, 0.0, 1.0).applyQuaternion(q);  // use it to figure out direction
        let vec2 = new THREE.Vector3(0.65, -0.76, 0.0).applyQuaternion(q); // ad hoc conversion between rpoly rotation and cadnano utils
        
        let new_position = move_along_vector(position , vec , n_bp)  // rpoly helix coordinates are defined in center of helix, cadnano utils have positions in the base of helix.
        
        for(const j of new_position.toArray())  {  // go through every coordinate to find the largest coordinate to figure out box size
            if (j > largest_size) {
                largest_size = j;
            }
        }

        // strand 0 is the scaffold and strand 1 is the staple
        let new_strands = generator.generate_or_sq(n_bp, undefined, new_position,vec, vec2);

        // cluster nucleotide by helix
        for(const i of [0,1])  {
            for(const nucleotide of new_strands[i]._nucleotides)  {
                nucleotide.cluster = n+1;
            }
        }

        // cut strand 1 into two equal lengh staple fragments for later connections
        let [fragment1, fragment2] = new_strands[1].cut_in_two(false);

        // store the fragments in this system for later connections
        staple_fragments.add_strand(fragment1);
        staple_fragments.add_strand(fragment2);

        scaffold_fragments.add_strand(new_strands[0]);
    }


    let output_system = new base.System(new THREE.Vector3(largest_size * 3.0, largest_size * 3.0, largest_size * 3.0));
    for(const n of rev_helix_connections)  {  // iterate through staple strand connections and connect the previously generated fragments
        let connect_from = n[0] * 2 - 1;
        let connect_to = n[1] * 2 - 2;
        let staple_strand = staple_fragments._strands[connect_from];
        staple_strand = staple_strand.append(staple_fragments._strands[connect_to]);
        output_system.add_strand(staple_strand)
    }

    let scaffold_strand = scaffold_fragments._strands[0];
    for(const n of fwd_helix_connections.slice(0,-1)) {
        let next_segment_adress = n[1]-1;
        let next_segment = scaffold_fragments._strands[next_segment_adress]
        scaffold_strand = scaffold_strand.append(next_segment);
    }

    scaffold_strand.make_circular();
    output_system.add_strand(scaffold_strand);

    return output_system;
}

export {loadRpoly}