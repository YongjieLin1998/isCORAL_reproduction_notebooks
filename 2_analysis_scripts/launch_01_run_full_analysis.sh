#!/bin/bash
# ==============================================================================
#        launch_01_run_full_analysis.sh
# ==============================================================================

# --- [1] Configuration Section ---
R_SCRIPT="2_analysis_scripts/01_run_full_analysis.R"
CONFIG_FILE="${1:-${ISCORAL_CONFIG_FILE:-config/<CONFIG_FILE>.yml}}"
R_RUNNER="${ISCORAL_R_RUNNER:-Rscript}"

# Automatically extract identifier from the config file for ledger records
CONFIG_FILENAME=$(basename "${CONFIG_FILE}")
SAMPLE_ID="${CONFIG_FILENAME%.*}" 
SH_SCRIPT=$(basename "$0")

# --- [2] Robust Working Directory Setup ---
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

if [ -z "${PROJECT_ROOT}" ] || [ ! -d "${PROJECT_ROOT}" ]; then
  echo "❌ CRITICAL ERROR: Failed to determine project root directory. Exiting."
  exit 1
fi
cd "${PROJECT_ROOT}"
echo ">>> Working directory set to: $(pwd)"

if [[ "${CONFIG_FILE}" == *"<CONFIG_FILE>"* ]] || [ ! -f "${CONFIG_FILE}" ]; then
  echo "CRITICAL ERROR: Provide an existing config file as the first argument or set ISCORAL_CONFIG_FILE."
  exit 2
fi

# --- [3] Archive & Ledger Initialization ---
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
ARCHIVE_ROOT="${PROJECT_ROOT}/8_run_archives/launch_01"
CURRENT_ARCHIVE="${ARCHIVE_ROOT}/run_${TIMESTAMP}_${SAMPLE_ID}"
MASTER_LOG="${PROJECT_ROOT}/8_run_archives/master_run_ledger.tsv"

echo ">>> [Archive] Initializing Run Archive & Ledger for ${SAMPLE_ID}..."
mkdir -p "${CURRENT_ARCHIVE}"

# Take a snapshot of the code and configuration
cp "$0" "${CURRENT_ARCHIVE}/launcher_snapshot.sh"
cp "${R_SCRIPT}" "${CURRENT_ARCHIVE}/r_script_snapshot.R"
if [ -f "${CONFIG_FILE}" ]; then
    cp "${CONFIG_FILE}" "${CURRENT_ARCHIVE}/config_snapshot.yml"
fi
env > "${CURRENT_ARCHIVE}/environment_vars.log"

if [ ! -f "${MASTER_LOG}" ]; then
    echo -e "Timestamp\tSample_ID\tLauncher_Script\tR_Script\tConfig_File\tArchive_Dir\tStatus" > "${MASTER_LOG}"
fi

echo -e "${TIMESTAMP}\t${SAMPLE_ID}\t${SH_SCRIPT}\t${R_SCRIPT}\t${CONFIG_FILE}\t${CURRENT_ARCHIVE}\tSTARTED" >> "${MASTER_LOG}"

# --- [4] Generate Background Wrapper Script (Robust Parameter Passing) ---
WRAPPER_SCRIPT="${CURRENT_ARCHIVE}/background_wrapper.sh"

# Note: Using quoted 'EOF' ensures that internal $ variables are not incorrectly expanded during generation
cat << 'EOF' > "${WRAPPER_SCRIPT}"
#!/bin/bash
# Enable job control, assigning an independent Process Group ID (PGID) to the core R script
set -m

# Receive safe parameters passed from the main script
W_PROJECT_ROOT="$1"
W_R_SCRIPT="$2"
W_CONFIG_FILE="$3"
W_CURRENT_ARCHIVE="$4"
W_TIMESTAMP="$5"
W_SAMPLE_ID="$6"
W_SH_SCRIPT="$7"
W_MASTER_LOG="$8"
W_R_RUNNER="$9"

# Define a function to update the ledger
update_ledger() {
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sed -i '' "s|${W_TIMESTAMP}\t${W_SAMPLE_ID}\t${W_SH_SCRIPT}\t${W_R_SCRIPT}\t${W_CONFIG_FILE}\t${W_CURRENT_ARCHIVE}\tSTARTED|${W_TIMESTAMP}\t${W_SAMPLE_ID}\t${W_SH_SCRIPT}\t${W_R_SCRIPT}\t${W_CONFIG_FILE}\t${W_CURRENT_ARCHIVE}\t${FINAL_STATUS}|g" "${W_MASTER_LOG}"
    else
        sed -i "s|${W_TIMESTAMP}\t${W_SAMPLE_ID}\t${W_SH_SCRIPT}\t${W_R_SCRIPT}\t${W_CONFIG_FILE}\t${W_CURRENT_ARCHIVE}\tSTARTED|${W_TIMESTAMP}\t${W_SAMPLE_ID}\t${W_SH_SCRIPT}\t${W_R_SCRIPT}\t${W_CONFIG_FILE}\t${W_CURRENT_ARCHIVE}\t${FINAL_STATUS}|g" "${W_MASTER_LOG}"
    fi
}

# 1. Launch core computation [Core modification: pass Archive directory as an argument to the R script]
"${W_R_RUNNER}" "${W_PROJECT_ROOT}/${W_R_SCRIPT}" "${W_PROJECT_ROOT}/${W_CONFIG_FILE}" "${W_CURRENT_ARCHIVE}" > "${W_CURRENT_ARCHIVE}/r_console_output.log" 2>&1 &
R_PID=$!

# 2. Process-group level signal interception (Trap)
trap 'echo ">>> User triggered kill signal. Terminating process tree..."; kill -TERM -$R_PID 2>/dev/null; FINAL_STATUS="KILLED_BY_USER"; update_ledger; exit 143' SIGTERM SIGINT

# 3. Wait for the process to terminate naturally
wait $R_PID
EXIT_CODE=$?

# 4. Record the final exit status
if [ $EXIT_CODE -eq 0 ]; then
    FINAL_STATUS="COMPLETED"
else
    FINAL_STATUS="FAILED_CODE_${EXIT_CODE}"
fi

update_ledger
EOF

chmod +x "${WRAPPER_SCRIPT}"


# --- [5] Submit the Job ---
echo ">>> [Compute] Preparing background job..."
echo "    R Script: ${R_SCRIPT}"
echo "    Config File: ${CONFIG_FILE}"
echo "    R Runner: ${R_RUNNER}"
echo "    Archive Dir: ${CURRENT_ARCHIVE}"

# Safely pass all necessary variables as arguments to the wrapper script
nohup bash "${WRAPPER_SCRIPT}" \
    "${PROJECT_ROOT}" \
    "${R_SCRIPT}" \
    "${CONFIG_FILE}" \
    "${CURRENT_ARCHIVE}" \
    "${TIMESTAMP}" \
    "${SAMPLE_ID}" \
    "${SH_SCRIPT}" \
    "${MASTER_LOG}" \
    "${R_RUNNER}" > /dev/null 2>&1 &

PID=$!
echo ""
echo "✅ Job successfully submitted to background!"
echo "📌 Process ID (PID) is: ${PID}"
echo "----------------------------------------------------------------------"
echo "👀 [Monitor Main Log]: tail -f ${CURRENT_ARCHIVE}/r_console_output.log"
echo "🔥 [Monitor Progress]: tail -f ${CURRENT_ARCHIVE}/realtime_progress.log"
echo "----------------------------------------------------------------------"
echo "🛑 Terminate this job securely with: kill ${PID}"

exit 0
