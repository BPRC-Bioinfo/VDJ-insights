document.addEventListener('DOMContentLoaded', function () {
    const simplifiedColumnsBtn = document.getElementById('applySimplifiedColumns');
    const allCheckboxes = document.querySelectorAll('#columnToggleForm input[name="columns"]');
    let isSimplified = true; // Track whether the button is in "Simplify" mode

    // Modal logic
    const modal = document.getElementById("columnToggleModal");
    const closeModalBtn = document.querySelector(".modal .close");
    const applyColumnsBtn = document.getElementById("applyColumns");

    document.querySelectorAll('.column-toggle-btn').forEach(btn => {
        btn.addEventListener('click', function () {
            modal.style.display = "flex"; // Ensure modal uses flex for centering
        });
    });

    closeModalBtn.onclick = function () {
        modal.style.display = "none";
    }

    window.onclick = function (event) {
        if (event.target == modal) {
            modal.style.display = "none";
        }
    }

    applyColumnsBtn.addEventListener('click', function () {
        applyColumnVisibility();
        modal.style.display = "none";
    });

    // Function to apply column visibility based on checked columns
    function applyColumnVisibility() {
        const checkboxes = document.querySelectorAll('#columnToggleForm input[name="columns"]:checked');
        const selectedColumns = Array.from(checkboxes).map(cb => cb.value);

        const tables = document.querySelectorAll('.report-table');
        tables.forEach(table => {
            const thElements = table.querySelectorAll('th');
            thElements.forEach(th => {
                const header = th.innerText.trim();
                const cellIndex = Array.from(th.parentElement.children).indexOf(th);

                // Always show the "Select" column
                if (header === "Select") {
                    th.style.display = ""; // Ensure the "Select" header is always shown
                    table.querySelectorAll('tr').forEach(tr => {
                        if (tr.children[cellIndex]) {
                            tr.children[cellIndex].style.display = ""; // Ensure the "Select" column is always shown
                        }
                    });
                } 
                // Handle all other columns including "Reference"
                else if (selectedColumns.includes(header)) {
                    th.style.display = ""; // Show the column if selected
                    table.querySelectorAll('tr').forEach(tr => {
                        if (tr.children[cellIndex]) {
                            tr.children[cellIndex].style.display = ""; // Show the column's cells
                        }
                    });
                } 
                // Hide columns not selected
                else {
                    th.style.display = "none"; // Hide the column
                    table.querySelectorAll('tr').forEach(tr => {
                        if (tr.children[cellIndex]) {
                            tr.children[cellIndex].style.display = "none"; // Hide the column's cells
                        }
                    });
                }
            });
        });
    }

    // Apply initial column visibility based on checked columns
    applyColumnVisibility();

    // Simplified columns button logic
    simplifiedColumnsBtn.addEventListener('click', function () {
        if (isSimplified) {
            // Simplify mode: check specific columns and switch to reset mode
            const columnsToCheck = ['Reference', 'Mismatches', 'Start coord', 'End coord', 'Region', 'Segment', 'Short name'];

            // First, uncheck all checkboxes
            allCheckboxes.forEach(checkbox => {
                checkbox.checked = false;
            });

            // Then, check only the specified columns
            columnsToCheck.forEach(columnName => {
                const checkbox = document.querySelector(`#columnToggleForm input[name="columns"][value="${columnName}"]`);
                if (checkbox) {
                    checkbox.checked = true;
                }
            });

            // Change button to "Reset" state
            simplifiedColumnsBtn.innerText = 'Reset';
            simplifiedColumnsBtn.classList.remove('simplified');
            simplifiedColumnsBtn.classList.add('reset');
            isSimplified = false; // Toggle to reset mode
        } else {
            // Reset mode: uncheck all checkboxes and switch back to simplify mode
            allCheckboxes.forEach(checkbox => {
                checkbox.checked = false;
            });

            // Change button back to "Simplify" state
            simplifiedColumnsBtn.innerText = 'Simplify';
            simplifiedColumnsBtn.classList.remove('reset');
            simplifiedColumnsBtn.classList.add('simplified');
            isSimplified = true; // Toggle back to simplify mode
        }
    });
    document.getElementById('go_button').addEventListener('click', function() {
        const sample_id = document.getElementById('sample_select').value;
        const annotation_type = document.getElementById('annotation_type').value;
        window.location.href = `/annotation/${annotation_type}/${sample_id}`;
    });    
});
