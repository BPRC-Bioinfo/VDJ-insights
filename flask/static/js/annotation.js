// annotation.js

document.addEventListener('DOMContentLoaded', function () {
    // Variables for simplified columns
    const simplifiedColumnsBtn = document.getElementById('applySimplifiedColumns');
    const allCheckboxes = document.querySelectorAll('#columnToggleForm input[name="columns"]');
    let isSimplified = true;

    // Modal logic for column toggle modal
    const columnModal = document.getElementById("columnToggleModal");
    const closeColumnModalBtn = columnModal.querySelector(".close");
    const applyColumnsBtn = document.getElementById("applyColumns");

    // Open column toggle modal
    document.querySelectorAll('.column-toggle-btn').forEach(btn => {
        btn.addEventListener('click', function () {
            columnModal.style.display = "flex";
        });
    });

    // Close column toggle modal when close button is clicked
    closeColumnModalBtn.onclick = function () {
        columnModal.style.display = "none";
    };

    // Apply column visibility when "Apply" button is clicked
    applyColumnsBtn.addEventListener('click', function () {
        applyColumnVisibility();
        columnModal.style.display = "none";
        updateRowCount(); // Update the row count after applying column visibility changes
    });

    // Modal logic for information modal
    const infoButton = document.getElementById('infoButton');
    const infoModal = document.getElementById('infoModal');
    const closeInfoModalBtn = infoModal.querySelector('.close');

    // Open info modal when info button is clicked
    if (infoButton) {
        infoButton.addEventListener('click', function () {
            infoModal.style.display = 'flex';
        });
    }

    // Close info modal when close button is clicked
    if (closeInfoModalBtn) {
        closeInfoModalBtn.addEventListener('click', function () {
            infoModal.style.display = 'none';
        });
    }

    // Close modals when clicking outside modal content
    window.addEventListener('click', function (event) {
        if (event.target === columnModal) {
            columnModal.style.display = "none";
        }
        if (event.target === infoModal) {
            infoModal.style.display = "none";
        }
    });

    // Close modals when pressing the Escape key
    window.addEventListener('keydown', function (event) {
        if (event.key === 'Escape') {
            if (columnModal.style.display === 'flex') {
                columnModal.style.display = 'none';
            }
            if (infoModal && infoModal.style.display === 'flex') {
                infoModal.style.display = 'none';
            }
        }
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
        const rowCountElement = document.getElementById('totalRowCount');
        if (rowCountElement) {
            rowCountElement.innerText = rowCount;
        }
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
        const annotationType = document.getElementById('annotation_type').value;
        let columnsToCheck = [];

        if (isSimplified) {
            if (annotationType === 'known') {
                columnsToCheck = ['Reference', 'Mismatches', 'Start coord', 'End coord', 'Region', 'Segment', 'Short name', 'Sample', 'Haplotype'];
            } else if (annotationType === 'novel') {
                columnsToCheck = ['Reference', 'Mismatches', 'Start coord', 'End coord', 'Region', 'Segment', 'Short name', 'Sample', 'Haplotype'];
            }

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

    // Event listener for the "Go" button
    document.getElementById('go_button').addEventListener('click', function () {
        const sample_id = document.getElementById('sample_select').value;
        const annotation_type = document.getElementById('annotation_type').value;
        window.location.href = `/annotation/${annotation_type}/${sample_id}`;
    });

    // Initial row count update on page load
    updateRowCount();

    // Hover effects for row count display
    const rowCountDisplay = document.getElementById('rowCountDisplay');

    // Apply hover styles to the element by default
    if (rowCountDisplay) {
        rowCountDisplay.classList.add('hover');

        rowCountDisplay.addEventListener('mouseenter', function () {
            const rowCount = document.getElementById('totalRowCount').innerText;
            rowCountDisplay.innerHTML = `<strong>Total rows:</strong> <span id="totalRowCount">${rowCount}</span>`;
        });

        rowCountDisplay.addEventListener('mouseleave', function () {
            const rowCount = document.getElementById('totalRowCount').innerText;
            rowCountDisplay.innerHTML = `<span id="totalRowCount">${rowCount}</span>`;
        });
    }
});
