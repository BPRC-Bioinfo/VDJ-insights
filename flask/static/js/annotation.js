document.addEventListener('DOMContentLoaded', function () {
    const simplifiedColumnsBtn = document.getElementById('applySimplifiedColumns');
    const allCheckboxes = document.querySelectorAll('#columnToggleForm input[name="columns"]');
    let isSimplified = true;

    // Modal logic
    const modal = document.getElementById("columnToggleModal");
    const closeModalBtn = document.querySelector(".modal .close");
    const applyColumnsBtn = document.getElementById("applyColumns");

    document.querySelectorAll('.column-toggle-btn').forEach(btn => {
        btn.addEventListener('click', function () {
            modal.style.display = "flex";
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
        updateRowCount(); // Update the row count after applying column visibility changes
    });

    // Function to update the total row count
    function updateRowCount() {
        const rows = document.querySelectorAll('#annotationTable tbody tr');
        let rowCount = 0;
        rows.forEach(row => {
            if (row.style.display !== 'none') {
                rowCount++;
            }
        });
        document.getElementById('totalRowCount').innerText = rowCount;
    }

    // Function to apply column visibility based on checked columns
    function applyColumnVisibility() {
        const checkboxes = document.querySelectorAll('#columnToggleForm input[name="columns"]');
        const selectedColumns = Array.from(checkboxes).filter(cb => cb.checked).map(cb => cb.value);

        const tables = document.querySelectorAll('.report-table');
        tables.forEach(table => {
            const thElements = table.querySelectorAll('th');
            thElements.forEach((th, index) => {
                const header = th.innerText.trim();

                if (header === "Select") {
                    th.style.display = "";
                    table.querySelectorAll('tr').forEach(tr => {
                        if (tr.children[index]) {
                            tr.children[index].style.display = "";
                        }
                    });
                } else if (selectedColumns.includes(header)) {
                    th.style.display = "";
                    table.querySelectorAll('tr').forEach(tr => {
                        if (tr.children[index]) {
                            tr.children[index].style.display = "";
                        }
                    });
                } else {
                    th.style.display = "none";
                    table.querySelectorAll('tr').forEach(tr => {
                        if (tr.children[index]) {
                            tr.children[index].style.display = "none";
                        }
                    });
                }
            });
        });
        updateRowCount(); // Update row count after adjusting column visibility
    }

    // Apply initial column visibility based on checked columns
    applyColumnVisibility();

    // Simplified columns button logic
    simplifiedColumnsBtn.addEventListener('click', function () {
        if (isSimplified) {
            const columnsToCheck = ['Reference', 'Mismatches', 'Start coord', 'End coord', 'Region', 'Segment', 'Short name', 'Sample', 'Haplotype'];

            allCheckboxes.forEach(checkbox => {
                checkbox.checked = columnsToCheck.includes(checkbox.value);
            });

            simplifiedColumnsBtn.innerText = 'Reset';
            simplifiedColumnsBtn.classList.remove('simplify');
            simplifiedColumnsBtn.classList.add('reset');
            isSimplified = false;
        } else {
            allCheckboxes.forEach(checkbox => {
                checkbox.checked = true;
            });

            simplifiedColumnsBtn.innerText = 'Simplify';
            simplifiedColumnsBtn.classList.remove('reset');
            simplifiedColumnsBtn.classList.add('simplify');
            isSimplified = true;
        }
        applyColumnVisibility(); // Apply changes and update row count
    });

    document.getElementById('go_button').addEventListener('click', function() {
        const sample_id = document.getElementById('sample_select').value;
        const annotation_type = document.getElementById('annotation_type').value;
        window.location.href = `/annotation/${annotation_type}/${sample_id}`;
    });

    // Initial row count update on page load
    updateRowCount();

    // Additional code to handle hover effects for row count display
    const rowCountDisplay = document.getElementById('rowCountDisplay');

    // Apply hover styles to the element by default
    rowCountDisplay.classList.add('hover');

    rowCountDisplay.addEventListener('mouseenter', function () {
        const rowCount = document.getElementById('totalRowCount').innerText;
        rowCountDisplay.innerHTML = `<strong>Total rows:</strong> <span id="totalRowCount">${rowCount}</span>`;
    });

    rowCountDisplay.addEventListener('mouseleave', function () {
        const rowCount = document.getElementById('totalRowCount').innerText;
        rowCountDisplay.innerHTML = `<span id="totalRowCount">${rowCount}</span>`;
    });
});
