var assert = require('assert');
//const THREE = require('three');
const tacoxdna = require('./dist/tacoxdna.js');

const fs = require('fs');

// Hide log output
tacoxdna.Logger.logFunction = (msg)=>{};

describe('Load cadnano', ()=>{
    const input = [fs.readFileSync('./tests/cadnano/input.json', 'utf8')];
    describe('Convert to oxDNA', ()=>{
        const [topology, conf] = tacoxdna.convertFromTo(
            input, 'cadnano', 'oxdna', {
                grid: 'sq',
                sequence: 'TATTCCCTCCCCCTACGATAAAGTGGTATTGTAGGGTCCAAGGATAAGTCTCGCACATAGCGAC',
        });
        const correctTop = fs.readFileSync('./tests/cadnano/correct_output.top', 'utf8');
        it('should should have correct topology', ()=>{
            assert.strictEqual(topology, correctTop);
        });
        const correctConf = fs.readFileSync('./tests/cadnano/correct_output.dat', 'utf8');
        it('should should have correct configuration', ()=>{
            assert.strictEqual(conf, correctConf);
        });
    });
    describe('Convert to oxView', ()=>{
        const output = tacoxdna.convertFromTo(
            input, 'cadnano', 'oxview', {
                grid: 'sq',
                sequence: 'TATTCCCTCCCCCTACGATAAAGTGGTATTGTAGGGTCCAAGGATAAGTCTCGCACATAGCGAC',
        });
        const correctOut = fs.readFileSync('./tests/cadnano/correct_output.oxview', 'utf8');
        it('should should have correct output', ()=>{
            assert.strictEqual(output, correctOut);
        });
    });
});

describe('Load PDB', ()=>{
    const input = [fs.readFileSync('./tests/pdb/input.pdb', 'utf8')];
    describe('Convert to oxDNA', ()=>{
        const [topology, conf] = tacoxdna.convertFromTo(
            input, 'pdb', 'oxdna', {strand_dir: '53'}
        );
        const correctTop = fs.readFileSync('./tests/pdb/correct_output.top', 'utf8');
        it('should should have correct topology', ()=>{
            assert.strictEqual(topology, correctTop);
        });
        const correctConf = fs.readFileSync('./tests/pdb/correct_output.dat', 'utf8');
        it('should should have correct configuration', ()=>{
            assert.strictEqual(conf, correctConf);
        });
    });
    describe('Convert to oxView', ()=>{
        const output = tacoxdna.convertFromTo(
            input, 'pdb', 'oxview', {strand_dir: '53'}
        );
        const correctOut = fs.readFileSync('./tests/pdb/correct_output.oxview', 'utf8');
        it('should should have correct output', ()=>{
            assert.strictEqual(output, correctOut);
        });
    });
});

describe('Load rpoly', ()=>{
    const input = [fs.readFileSync('./tests/rpoly/input.rpoly', 'utf8')];
    const seq = 'TACAATCCGTACGACGAAACAAGTTAAATAAGATAAACAATGTTGTTTCATCCCACGTAGTAGTTAAACACGTTTGGCAGCCGCCCTGCTAGCCCCCTTATTTCGACGTCGATGTCGCAACTGAATCTCCATGCCAGCTGTTACGGGTGAGGTTAGCCACAGTCAGATGGATATATCAGGAGAATCTGCCTGAGTCCCTCCGGTCTACAAGGTCTGAAAAAATATAGGGTCCAAGGATAAGTCTCGCACATAGCGACAGACGCATTTTCAGAACAACCGCATATTCCAATGTTATGGTGAAATAGCATCCCCCTCCCTTATACCAATATTTTAGCCG';

    describe('Convert to oxDNA', ()=>{
        const [topology, conf] = tacoxdna.convertFromTo(
            input, 'rpoly', 'oxdna', {
                sequence: seq
            }
        );

        const correctTop = fs.readFileSync('./tests/rpoly/correct_output.top', 'utf8');
        it('should should have correct topology', ()=>{
            assert.strictEqual(topology, correctTop);
        });

        const correctConf = fs.readFileSync('./tests/rpoly/correct_output.dat', 'utf8');
        it('should should have correct configuration', ()=>{
            assert.strictEqual(conf, correctConf);
        });

    });

    describe('Convert to oxView', ()=>{
        const output = tacoxdna.convertFromTo(
            input, 'rpoly', 'oxview', {
                sequence: seq
            }
        );
        const correctOut = fs.readFileSync('./tests/rpoly/correct_output.oxview', 'utf8');
        it('should should have correct output', ()=>{
            assert.strictEqual(output, correctOut);
        });
    });
});

describe('Load tiamat', ()=>{
    const input = [fs.readFileSync('./tests/tiamat/input.dnajson', 'utf8')];

    describe('Convert to oxDNA', ()=>{
        const [topology, conf] = tacoxdna.convertFromTo(
            input, 'tiamat', 'oxdna', {
                tiamat_version: 2,
                isDNA: true,
                default_val: "R"
            }
        );

        const correctTop = fs.readFileSync('./tests/tiamat/correct_output.top', 'utf8');
        it('should should have correct topology', ()=>{
            assert.strictEqual(topology, correctTop);
        });

        const correctConf = fs.readFileSync('./tests/tiamat/correct_output.dat', 'utf8');
        it('should should have correct configuration', ()=>{
            assert.strictEqual(conf, correctConf);
        });

    });

    describe('Convert to oxView', ()=>{
        const output = tacoxdna.convertFromTo(
            input, 'tiamat', 'oxview', {
                tiamat_version: 2,
                isDNA: true,
                default_val: "R"
            }
        );
        const correctOut = fs.readFileSync('./tests/tiamat/correct_output.oxview', 'utf8');
        it('should should have correct output', ()=>{
            assert.strictEqual(output, correctOut);
        });
    });
});
