var assert = require('assert');
//const THREE = require('three');
const tacoxdna = require('./dist/tacoxdna.js');

const fs = require('fs');

describe('Load cadnano', function() {
    const cadnanoFile = fs.readFileSync('./tests/cadnano/init.json').toString();
    describe('Convert to oxDNA', function() {
        const [topology, conf] = tacoxdna.convertFromTo(
            [cadnanoFile], 'cadnano', 'oxdna', {
                grid: 'sq',
                sequence: 'TATTCCCTCCCCCTACGATAAAGTGGTATTGTAGGGTCCAAGGATAAGTCTCGCACATAGCGAC',
        });
        const correctTop = fs.readFileSync('./tests/cadnano/correct_output.top').toString();
        it('should should have correct topology', function() {
            assert.strictEqual(topology, correctTop);
        });
        const correctConf = fs.readFileSync('./tests/cadnano/correct_output.dat').toString();
        it('should should have correct configuration', function() {
            assert.strictEqual(conf, correctConf);
        });
    });
    describe('Convert to oxView', function() {
        const output = tacoxdna.convertFromTo(
            [cadnanoFile], 'cadnano', 'oxview', {
                grid: 'sq',
                sequence: 'TATTCCCTCCCCCTACGATAAAGTGGTATTGTAGGGTCCAAGGATAAGTCTCGCACATAGCGAC',
        });
        const correctOut = fs.readFileSync('./tests/cadnano/correct_output.oxview').toString();
        it('should should have correct output', function() {
            assert.strictEqual(output, correctOut);
        });
    });
});
