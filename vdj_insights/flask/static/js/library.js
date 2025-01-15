document.addEventListener('DOMContentLoaded', function () {
    // Select All functionality for checkboxes in the "Select" column
    const selectAllBtn = document.getElementById('selectAllBtn');
    let isSelectAll = true;

    if (selectAllBtn) {  // Only add event listener if button exists
        selectAllBtn.addEventListener('click', function () {
            // Select only checkboxes in rows that are visible
            const visibleRows = document.querySelectorAll('.report-table tbody tr');
            visibleRows.forEach(row => {
                if (row.style.display !== 'none') {
                    const checkbox = row.querySelector('.row-checkbox');
                    if (checkbox) {
                        checkbox.checked = isSelectAll;
                    }
                }
            });

            // Toggle the button text between 'Select All' and 'Deselect All'
            selectAllBtn.textContent = isSelectAll ? 'Deselect All' : 'Select All';
            isSelectAll = !isSelectAll;
        });
    }

    // Save to Library functionality
    const saveToLibraryBtn = document.getElementById('saveToLibraryBtn');
    if (saveToLibraryBtn) {  // Only add event listener if button exists
        saveToLibraryBtn.addEventListener('click', function () {
            const selectedRows = document.querySelectorAll('.row-checkbox:checked');

            if (selectedRows.length === 0) {
                alert('No rows selected');
                return;
            }

            let newSequences = {};
            const table = document.querySelector('.report-table');
            const headers = Array.from(table.querySelectorAll('th'));

            const oldNameLikeIndex = headers.findIndex(header => header.innerText.trim() === "Old name-like");
            const startCoordIndex = headers.findIndex(header => header.innerText.trim() === "Start coord");
            const endCoordIndex = headers.findIndex(header => header.innerText.trim() === "End coord");
            const sampleIndex = headers.findIndex(header => header.innerText.trim() === "Sample");
            const sequenceIndex = headers.findIndex(header => header.innerText.trim() === "Old name-like seq");

            selectedRows.forEach(checkbox => {
                const row = checkbox.closest('tr');
                const rowCells = Array.from(row.children);

                const oldNameLike = rowCells[oldNameLikeIndex].innerText;
                const startCoord = rowCells[startCoordIndex].innerText;
                const endCoord = rowCells[endCoordIndex].innerText;
                const sample = rowCells[sampleIndex].innerText;
                const sequence = rowCells[sequenceIndex].innerText;

                const fastaHeader = `>${oldNameLike}_${startCoord}_${endCoord}_${sample}`;
                newSequences[fastaHeader] = sequence;
            });

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
    }

    // Remove from Library functionality
    const removeFromLibraryBtn = document.getElementById('removeFromLibraryBtn');
    if (removeFromLibraryBtn) {  // Only add event listener if button exists
        removeFromLibraryBtn.addEventListener('click', function () {
            const selectedRows = document.querySelectorAll('.row-checkbox:checked');

            if (selectedRows.length === 0) {
                alert('No rows selected');
                return;
            }

            let sequencesToRemove = [];
            const table = document.querySelector('.report-table');
            const headers = Array.from(table.querySelectorAll('th'));
            const headerIndex = headers.findIndex(header => header.innerText.trim() === "Header");

            selectedRows.forEach(checkbox => {
                const row = checkbox.closest('tr');
                const rowCells = Array.from(row.children);
                const sequenceHeader = rowCells[headerIndex].innerText;
                sequencesToRemove.push(sequenceHeader);
            });

            fetch('/remove_sequences', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({ sequences: sequencesToRemove })
            })
            .then(response => response.json())
            .then(data => {
                if (data.message) {
                    alert(data.message);
                    // Remove the selected rows from the table
                    selectedRows.forEach(checkbox => {
                        const row = checkbox.closest('tr');
                        row.remove();
                    });
                } else {
                    alert('Error removing sequences');
                }
            })
            .catch(error => console.error('Error:', error));
        });
    }
    const downloadBtn = document.getElementById('downloadFastaBtn');

    if (downloadBtn) {
        downloadBtn.addEventListener('click', function (event) {
            event.preventDefault(); // Prevent any default button action

            fetch('/download_fasta')
                .then(response => {
                    if (!response.ok) {
                        throw new Error('Failed to download the file.');
                    }
                    return response.blob();
                })
                .then(blob => {
                    const url = window.URL.createObjectURL(blob);
                    const a = document.createElement('a');
                    a.style.display = 'none';
                    a.href = url;
                    a.download = 'library.fasta'; // The default filename for the download
                    document.body.appendChild(a);
                    a.click();
                    window.URL.revokeObjectURL(url); // Clean up
                })
                .catch(error => {
                    console.error('Error downloading FASTA:', error);
                    alert('An error occurred while downloading the file.');
                });
        });
    }
});
