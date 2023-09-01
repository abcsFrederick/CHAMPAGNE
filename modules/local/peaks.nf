
process SICER {
    label 'peaks'

    input:
        tuple val(meta), path(chip), path(input)

}
