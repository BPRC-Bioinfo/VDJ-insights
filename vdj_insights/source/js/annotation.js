document.addEventListener('DOMContentLoaded', function () {
    // Toggle additional content sections
    const showMoreButtons = document.querySelectorAll('.modal-btn');
    showMoreButtons.forEach(button => {
        button.addEventListener('click', function () {
            const section = document.getElementById(this.dataset.section);
            if (section) {
                section.classList.toggle('active');
                this.textContent = section.classList.contains('active') ? 'Show Less' : 'Show All';
            }
        });
    });

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

    // Select All functionality for checkboxes in the "Select" column
    const selectAllBtn = document.getElementById('selectAllBtn');
    let isSelectAll = true; // Track state of select/deselect all

    selectAllBtn.addEventListener('click', function () {
        const checkboxes = document.querySelectorAll('.row-checkbox');
        checkboxes.forEach(checkbox => {
            checkbox.checked = isSelectAll; // Check or uncheck all checkboxes
        });
        // Toggle the button text and the select/deselect state
        selectAllBtn.textContent = isSelectAll ? 'Deselect All' : 'Select All';
        isSelectAll = !isSelectAll; // Toggle the select/deselect state
    });

    // Save to Library functionality
    const saveToLibraryBtn = document.getElementById('saveToLibraryBtn');
    saveToLibraryBtn.addEventListener('click', function () {
        const selectedRows = document.querySelectorAll('.row-checkbox:checked');

        if (selectedRows.length === 0) {
            alert('No rows selected');
            return;
        }

        let newSequences = {};

        // Get column indexes by their header names
        const table = document.querySelector('.report-table');
        const headers = Array.from(table.querySelectorAll('th'));

        const oldNameLikeIndex = headers.findIndex(header => header.innerText.trim() === "Old name-like");
        const startCoordIndex = headers.findIndex(header => header.innerText.trim() === "Start coord");
        const endCoordIndex = headers.findIndex(header => header.innerText.trim() === "End coord");
        const sampleIndex = headers.findIndex(header => header.innerText.trim() === "Sample");
        const sequenceIndex = headers.findIndex(header => header.innerText.trim() === "Old name-like seq");

        // Loop through each selected row
        selectedRows.forEach(checkbox => {
            const row = checkbox.closest('tr'); // Get the row where the checkbox is checked

            const rowCells = Array.from(row.children);

            // Extract data from the correct columns using their indexes
            const oldNameLike = rowCells[oldNameLikeIndex].innerText;
            const startCoord = rowCells[startCoordIndex].innerText;
            const endCoord = rowCells[endCoordIndex].innerText;
            const sample = rowCells[sampleIndex].innerText;
            const sequence = rowCells[sequenceIndex].innerText;

            const fastaHeader = `>${oldNameLike}_${startCoord}_${endCoord}_${sample}`;

            // Add to newSequences object
            newSequences[fastaHeader] = sequence;
        });

        // Send the selected sequences to the backend
        fetch('/save_sequences', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ sequences: newSequences })
        })
        .then(response => response.json())
        .then(data => {
            if (data.message) {
                alert(data.message);
            } else {
                alert('Error saving sequences');
            }
        })
        .catch(error => console.error('Error:', error));
    });
});
