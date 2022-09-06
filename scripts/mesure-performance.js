
const { promisify } = require('util')
const { parse } = require("csv-parse");
const fs = require('fs');

const exec = promisify(require('child_process').exec)

const csvMock = ``

const ITERATIONS_LIMIT = 10

const LIKWID_COMMAND = 'likwid-perfctr'
const FIRST_FLAGS = '-O -C 1 -g'
const SECOND_FLAGS = '-m'
const PROGRAM = 'invmat'

const L3 = 'L3'
const L2CACHE = 'L2CACHE'
const FLOPS_DP = 'FLOPS_DP'
const AVX_FLOPS_DP = 'FLOPS_DP'

const groups = [L3, L2CACHE, FLOPS_DP, AVX_FLOPS_DP]
const sizes = [64, 100, 128, 1024, 2000, 2048]

const getFileContents = async (filepath) => {
    const data = [];
    return new Promise(function(resolve, reject) {
        fs.createReadStream(filepath)
            .pipe(parse({delimiter: ',', relax_column_count: true}))
            .on('error', error => reject(error))
            .on('data', row => data.push(row))
            .on('end', () => {
                resolve(data);
            });
    });
}

const parseKey = async (key) => {
    //const data = await getFileContents('test.csv')
    const data = await getFileContents('output.csv')

    let linearSystemCalculationL3 = 0;
    let i;
    for (i = 0; i < data.length; i++) {
        if (data[i][0] === key) {
            linearSystemCalculationL3 = +(data[i][1])
            break;
        }
    }

    let residueL3 = 0;
    for (; i < data.length; i++) {
        if (data[i][0] === key) {
            residueL3 = +(data[i][1])
        }
    }

    return {
        'linearSystem': linearSystemCalculationL3,
        'residue': residueL3
    }
}

const parseL3 = async () => {
    return parseKey('L3 bandwidth [MBytes/s]')
}

const parseL2Cache = () => {
    return parseKey('L2 miss ratio')
}

const parseFlopsDP = () => {
    return parseKey('DP MFLOP/s')
}

const parseAVXFlopsDP = () => {
    return parseKey('AVX DP MFLOP/s')
}

const parsers = {
    L3: parseL3,
    L2CACHE: parseL2Cache,
    FLOPS_DP: parseFlopsDP,
    AVX_FLOPS_DP: parseAVXFlopsDP
}

const buildCommand = (group, size) => {
    return `${LIKWID_COMMAND} ${FIRST_FLAGS} ${group} ${SECOND_FLAGS} ./${PROGRAM} -r ${size} -i ${ITERATIONS_LIMIT} -s output.csv`
}

const execMock = async (command) => {
    console.log('Executing command: ' + command)
    return {
        stdout: csvMock
    }
}

const run = async (group, size, parser, mockExecution) => {
    const command = buildCommand(group, size)
    try {
        if (mockExecution) {
            await execMock(command)
            return parser()
        }
        console.log('Executing command: ' + command)
        await exec(command)
        return parser()
    } catch (e) {
        console.error(e)
        return parser()
    }
}

function Result(size, indicators) {
    this.matrixSize = size
    for (const indicator of indicators) {
        this[indicator.group] = indicator.value
    }
}

const main = async () => {

    const mockExecution = process.argv[2] === 'mock'
    const fixedSize = parseInt(process.argv[3])
    let actualSizes = sizes;
    if (fixedSize > 0) {
        actualSizes = [fixedSize]
    }

    const results = await Promise.all(actualSizes.map(async (size) => {

        const indicators = await Promise.all(groups.map(async group => ({
            size: size,
            group: group,
            values: await run(group, size, parsers[group], mockExecution)
        })))

        const linearSystemIndicators = indicators.map(indicator => {
            return {
                size: indicator.size,
                group: indicator.group,
                value: indicator.values.linearSystem
            }
        })

        const residueIndicators = indicators.map(indicator => {
            return {
                size: indicator.size,
                group: indicator.group,
                value: indicator.values.residue
            }
        })

        return {
            linearSystemResult: new Result(size, linearSystemIndicators),
            residueResult: new Result(size, residueIndicators),
        }
    }))

    const linearSystemResults = results.map(result => result.linearSystemResult)
    const residueResults = results.map(result => result.residueResult)

    console.log('\nLINEAR SYSTEM RESULTS:\n')

    console.table(linearSystemResults)

    console.log('\nRESIDUE RESULTS:\n')

    console.table(residueResults)
}

main()