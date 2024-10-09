// sample_details.js

document.addEventListener('DOMContentLoaded', function () {
    // Get the sample ID from the summary header
    const sample_id = document.querySelector('.sample-summary-header h1')?.textContent;

    // Get all collapsible elements
    const collapsibles = document.querySelectorAll('.collapsible');

    // Add click event listener to each collapsible
    collapsibles.forEach(function (collapsible) {
        collapsible.addEventListener('click', function () {
            // Toggle the current section
            this.classList.toggle('active');
            const content = this.nextElementSibling;

            if (this.classList.contains('active')) {
                content.style.display = 'block';
            } else {
                content.style.display = 'none';
            }
        });
    });

    // Download Region Data as CSV
    const downloadRegionButton = document.getElementById('downloadRegionData');
    if (downloadRegionButton) {
        downloadRegionButton.addEventListener('click', function () {
            downloadTableAsCSV('regionDataTable', `${sample_id}_region_data.csv`);
        });
    }

    // Download Additional Data as CSV
    const downloadAdditionalButton = document.getElementById('downloadAdditionalData');
    if (downloadAdditionalButton) {
        downloadAdditionalButton.addEventListener('click', function () {
            downloadTableAsCSV('additionalDataTable', `${sample_id}_additional_data.csv`);
        });
    }

    // Function to download table data as CSV
    function downloadTableAsCSV(tableId, filename) {
        const table = document.getElementById(tableId);
        if (!table) return;

        const rows = table.querySelectorAll('tr');
        let csvContent = '';

        rows.forEach(function (row) {
            const cols = row.querySelectorAll('th, td');
            const rowData = [];
            cols.forEach(function (col) {
                rowData.push('"' + col.innerText.replace(/"/g, '""') + '"');
            });
            csvContent += rowData.join(',') + '\n';
        });

        // Create a Blob and trigger download
        const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
        const url = URL.createObjectURL(blob);

        const link = document.createElement('a');
        link.href = url;
        link.download = filename;
        link.style.display = 'none';
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    }
});
