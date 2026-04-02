#!/bin/bash
#SBATCH --job-name=build_container
#SBATCH --output=logs/build_container_%j.out
#SBATCH --error=logs/build_container_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal

################################################################################
# Apptainer Container Build Script
# Builds coloc-pipeline.sif from coloc-pipeline.def
################################################################################

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Use local scratch for temporary files (faster I/O during build)
export APPTAINER_TMPDIR=/localscratch/${USER}
mkdir -p $APPTAINER_TMPDIR

# Clean up local scratch on exit
trap "log_info 'Cleaning up temporary files...'; rm -rf ${APPTAINER_TMPDIR}/sbuild-* ${APPTAINER_TMPDIR}/build-temp-* 2>/dev/null || true" EXIT

# Configuration
PIPELINE_DIR="/home/itg/oleg.vlasovets/projects/coloc-pipeline"
DEF_FILE="${PIPELINE_DIR}/coloc-pipeline.def"
SIF_FILE="${PIPELINE_DIR}/coloc-pipeline.sif"

# Functions
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

check_prerequisites() {
    log_info "Checking prerequisites..."
    
    # Check if definition file exists
    if [[ ! -f "${DEF_FILE}" ]]; then
        log_error "Definition file not found: ${DEF_FILE}"
        log_error "Please create coloc-pipeline.def before running this script"
        exit 1
    fi
    
    # Check if apptainer is available
    if ! command -v apptainer &> /dev/null; then
        log_error "Apptainer not found in PATH"
        log_error "Please load the apptainer module or ensure it's installed"
        exit 1
    fi
    
    # Check apptainer version
    APPTAINER_VERSION=$(apptainer --version 2>&1 || echo "unknown")
    log_info "Apptainer version: ${APPTAINER_VERSION}"
    
    # Check if user has fakeroot capability
    if apptainer version 2>&1 | grep -q "fakeroot"; then
        log_info "Fakeroot capability: available"
    else
        log_info "Fakeroot capability: checking..."
    fi
    
    log_info "Prerequisites OK"
}

remove_existing_image() {
    if [[ -f "${SIF_FILE}" ]]; then
        log_info "Removing existing image: ${SIF_FILE}"
        EXISTING_SIZE=$(du -sh "${SIF_FILE}" | cut -f1)
        log_info "Previous image size: ${EXISTING_SIZE}"
        rm -f "${SIF_FILE}"
    fi
}

build_container() {
    log_info "=========================================="
    log_info "Building Apptainer Container"
    log_info "=========================================="
    log_info "Definition file: ${DEF_FILE}"
    log_info "Output image: ${SIF_FILE}"
    log_info "Temporary directory: ${APPTAINER_TMPDIR}"
    log_info "CPUs available: ${SLURM_CPUS_PER_TASK:-4}"
    log_info "Memory available: ${SLURM_MEM_PER_NODE:-16384}MB"
    log_info ""
    
    # Change to pipeline directory
    cd "${PIPELINE_DIR}" || exit 1
    
    # Build start time
    BUILD_START=$(date +%s)
    
    # Build the container with fakeroot
    log_info "Starting build (this may take 30-60 minutes)..."
    log_info "Command: apptainer build --fakeroot --ignore-fakeroot-command ${SIF_FILE} ${DEF_FILE}"
    log_info ""
    
    if apptainer build --fakeroot --ignore-fakeroot-command "${SIF_FILE}" "${DEF_FILE}"; then
        # Build end time
        BUILD_END=$(date +%s)
        BUILD_TIME=$((BUILD_END - BUILD_START))
        BUILD_MINUTES=$((BUILD_TIME / 60))
        BUILD_SECONDS=$((BUILD_TIME % 60))
        
        log_info "=========================================="
        log_info "Build completed successfully!"
        log_info "Build time: ${BUILD_MINUTES}m ${BUILD_SECONDS}s"
        log_info "=========================================="
        return 0
    else
        BUILD_EXIT_CODE=$?
        log_error "=========================================="
        log_error "Build failed with exit code ${BUILD_EXIT_CODE}"
        log_error "=========================================="
        log_error "Check the log files:"
        log_error "  STDOUT: logs/build_container_${SLURM_JOB_ID}.out"
        log_error "  STDERR: logs/build_container_${SLURM_JOB_ID}.err"
        return "${BUILD_EXIT_CODE}"
    fi
}

test_container() {
    log_info ""
    log_info "=========================================="
    log_info "Testing Container"
    log_info "=========================================="
    
    # Run %test block
    log_info "Running built-in tests (%%test block)..."
    if apptainer test "${SIF_FILE}"; then
        log_info "Built-in tests PASSED"
    else
        log_error "Built-in tests FAILED"
        return 1
    fi
    
    # Quick smoke tests
    log_info ""
    log_info "Running smoke tests..."
    
    # Test 1: Check if R is available
    log_info "  Test 1: R availability"
    if apptainer exec "${SIF_FILE}" R --version > /dev/null 2>&1; then
        R_VERSION=$(apptainer exec "${SIF_FILE}" R --version 2>&1 | head -1)
        log_info "    ✓ R found: ${R_VERSION}"
    else
        log_error "    ✗ R not found or not working"
        return 1
    fi
    
    # Test 2: Check if Snakemake is available
    log_info "  Test 2: Snakemake availability"
    if apptainer exec "${SIF_FILE}" snakemake --version > /dev/null 2>&1; then
        SNAKEMAKE_VERSION=$(apptainer exec "${SIF_FILE}" snakemake --version 2>&1)
        log_info "    ✓ Snakemake found: ${SNAKEMAKE_VERSION}"
    else
        log_error "    ✗ Snakemake not found or not working"
        return 1
    fi
    
    # Test 3: Check if bcftools is available
    log_info "  Test 3: bcftools availability"
    if apptainer exec "${SIF_FILE}" bcftools --version > /dev/null 2>&1; then
        BCFTOOLS_VERSION=$(apptainer exec "${SIF_FILE}" bcftools --version 2>&1 | head -1)
        log_info "    ✓ bcftools found: ${BCFTOOLS_VERSION}"
    else
        log_error "    ✗ bcftools not found or not working"
        return 1
    fi
    
    log_info ""
    log_info "All smoke tests PASSED"
    log_info "=========================================="
    return 0
}

show_summary() {
    log_info ""
    log_info "=========================================="
    log_info "Build Summary"
    log_info "=========================================="
    
    # Image size
    IMAGE_SIZE=$(du -sh "${SIF_FILE}" | cut -f1)
    log_info "Image file: ${SIF_FILE}"
    log_info "Image size: ${IMAGE_SIZE}"
    
    # Image info
    log_info ""
    log_info "Container details:"
    apptainer inspect "${SIF_FILE}" 2>/dev/null | grep -E "^(Author|Version|Description):" || true
    
    log_info ""
    log_info "=========================================="
    log_info "Next steps:"
    log_info "=========================================="
    log_info "1. Test interactively:"
    log_info "   apptainer shell ${SIF_FILE}"
    log_info ""
    log_info "2. Run a command:"
    log_info "   apptainer exec ${SIF_FILE} R --version"
    log_info ""
    log_info "3. Submit a job using the container:"
    log_info "   sbatch scripts/stage2_overlaps_container.sh"
    log_info ""
    log_info "See CONTAINER.md for full usage documentation"
    log_info "=========================================="
}

# Main execution
main() {
    log_info "=========================================="
    log_info "Apptainer Container Build"
    log_info "=========================================="
    log_info "SLURM Job ID: ${SLURM_JOB_ID:-N/A}"
    log_info "Node: ${SLURM_NODELIST:-$(hostname)}"
    log_info "=========================================="
    log_info ""
    
    # Check prerequisites
    check_prerequisites
    
    # Remove existing image
    remove_existing_image
    
    # Build the container
    if ! build_container; then
        EXIT_CODE=$?
        log_error ""
        log_error "Build failed. See logs for details:"
        log_error "  ${PIPELINE_DIR}/logs/build_container_${SLURM_JOB_ID}.out"
        log_error "  ${PIPELINE_DIR}/logs/build_container_${SLURM_JOB_ID}.err"
        exit "${EXIT_CODE}"
    fi
    
    # Test the container
    if ! test_container; then
        log_error ""
        log_error "Container tests failed!"
        log_error "Image was built but may not be functional"
        exit 1
    fi
    
    # Show summary
    show_summary
    
    log_info ""
    log_info "=========================================="
    log_info "BUILD SUCCESSFUL!"
    log_info "=========================================="
    exit 0
}

# Run main function
main "$@"
