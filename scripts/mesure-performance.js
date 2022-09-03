
const { promisify } = require('util')

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

const groups = [L3, L2CACHE, FLOPS_DP]
const sizes = [64, 100, 128, 1024, 2000, 2048]

const parseL3 = (csv) => {
    return 123.32
}

const parseL2Cache = (csv) => {
    return 123.32
}

const parseFlopsDP = (csv) => {
    return 123.32
}

const parsers = {
    L3: parseL3,
    L2CACHE: parseL2Cache,
    FLOPS_DP: parseFlopsDP,
}

const buildCommand = (group, size) => {
    return `${LIKWID_COMMAND} ${FIRST_FLAGS} ${group} ${SECOND_FLAGS} ./${PROGRAM} -r ${size} -i ${ITERATIONS_LIMIT}`
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
            const resultMock = await execMock(command)
            return parser(resultMock.stdout)
        }
        const result = await exec(command)
        return parser(result.stdout)
    } catch (e) {
        console.error(e)
        return parser(e.stdout)
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

    const results = await Promise.all(sizes.map(async (size) => {

        const indicators = await Promise.all(groups.map(async group => ({
            group: group,
            value: await run(group, size, parsers[group], mockExecution)
        })))

        return new Result(size, indicators)
    }))

    console.log('\nRESULTS:\n')

    console.table(results)
}

main()