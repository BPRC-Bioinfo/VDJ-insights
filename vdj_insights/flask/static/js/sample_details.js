document.addEventListener('DOMContentLoaded', function () {
    // Get the sample ID from the summary header
    const sample_id = document.querySelector('.sample-summary-header h1')?.textContent;

    // Initialize collapsible sections
    const collapsibles = document.querySelectorAll('.collapsible');
    collapsibles.forEach(function (collapsible) {
        collapsible.classList.add('active');
        const content = collapsible.nextElementSibling;
        content.style.display = 'block';

        collapsible.addEventListener('click', function () {
            this.classList.toggle('active');
            content.style.display = this.classList.contains('active') ? 'block' : 'none';
        });
    });

    // Download Region Data
    const downloadRegionButton = document.getElementById('downloadRegionData');
    if (downloadRegionButton) {
        downloadRegionButton.addEventListener('click', function () {
            const regionTable = document.getElementById('regionDataTable');
            if (!regionTable) {
                alert('No region data available for download.');
                return;
            }
            downloadTableAsCSV(regionTable, `${sample_id}_region_data.csv`);
        });
    }

    // Function to download a table as CSV
    function downloadTableAsCSV(table, filename) {
        let csvContent = '';
        const rows = table.querySelectorAll('tr');
        rows.forEach(function (row) {
            const cols = row.querySelectorAll('th, td');
            const rowData = Array.from(cols).map(col => {
                let data = col.innerText.replace(/"/g, '""');
                if (/[,"\n]/.test(data)) {
                    data = `"${data}"`;
                }
                return data;
            });
            csvContent += rowData.join(',') + '\n';
        });

        // Trigger download
        const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
        const link = document.createElement('a');
        link.href = URL.createObjectURL(blob);
        link.download = filename;
        link.click();
    }

    // Download Annotation Data as Excel
    const downloadAnnotationButton = document.getElementById('downloadAnnotationData');
    if (downloadAnnotationButton) {
        downloadAnnotationButton.addEventListener('click', function () {
            const annotationSection = document.querySelector('.statistics-container');
            if (!annotationSection) {
                alert('No annotation data available for download.');
                return;
            }
            downloadAnnotationSectionAsExcel(annotationSection, `${sample_id}_annotation_data.xlsx`);
        });
    }

    // Function to collect known and novel annotation data and download it as Excel
    function downloadAnnotationSectionAsExcel(section, filename) {
        const workbook = XLSX.utils.book_new();
        
        // Prepare data for Known Annotations sheet
        const knownAnnotationsData = [['Known Annotations']];
        const knownSummary = collectAnnotationSummary(section, 'Known Annotations');
        if (knownSummary) {
            knownAnnotationsData.push(...knownSummary);
            const knownWorksheet = XLSX.utils.aoa_to_sheet(knownAnnotationsData);
            XLSX.utils.book_append_sheet(workbook, knownWorksheet, 'Known Annotations');
        }

        // Prepare data for Novel Annotations sheet
        const novelAnnotationsData = [['Novel Annotations']];
        const novelSummary = collectAnnotationSummary(section, 'Novel Annotations');
        if (novelSummary) {
            novelAnnotationsData.push(...novelSummary);
            const novelWorksheet = XLSX.utils.aoa_to_sheet(novelAnnotationsData);
            XLSX.utils.book_append_sheet(workbook, novelWorksheet, 'Novel Annotations');
        }

        // Trigger Excel file download
        XLSX.writeFile(workbook, filename);
    }

    // Function to collect annotation summary data
    function collectAnnotationSummary(section, headerText) {
        const headers = section.querySelectorAll('h3');
        for (let header of headers) {
            if (header.textContent.trim() === headerText) {
                const summaryData = [];
                const summaryTable = header.nextElementSibling; // Total Annotations, etc.

                // Collect summary stats
                const rows = summaryTable.querySelectorAll('tr');
                rows.forEach(function (row) {
                    const cols = row.querySelectorAll('th, td');
                    const rowData = Array.from(cols).map(col => col.innerText);
                    summaryData.push(rowData);
                });

                // Find and add "Segments Distribution"
                const segmentsHeader = findHeader(section, 'Segments Distribution');
                if (segmentsHeader) {
                    summaryData.push(['Segments Distribution']);
                    const segmentsTable = segmentsHeader.nextElementSibling;
                    const segmentRows = segmentsTable.querySelectorAll('tr');
                    segmentRows.forEach(function (row) {
                        const cols = row.querySelectorAll('th, td');
                        const rowData = Array.from(cols).map(col => col.innerText);
                        summaryData.push(rowData);
                    });
                }

                // Find and add "Functions Distribution"
                const functionsHeader = findHeader(section, 'Functions Distribution');
                if (functionsHeader) {
                    summaryData.push(['Functions Distribution']);
                    const functionsTable = functionsHeader.nextElementSibling;
                    const functionRows = functionsTable.querySelectorAll('tr');
                    functionRows.forEach(function (row) {
                        const cols = row.querySelectorAll('th, td');
                        const rowData = Array.from(cols).map(col => col.innerText);
                        summaryData.push(rowData);
                    });
                }

                return summaryData;
            }
        }
        return null;
    }

    // Function to find headers like "Segments Distribution" or "Functions Distribution"
    function findHeader(section, headerText) {
        const headers = section.querySelectorAll('h4');
        for (let header of headers) {
            if (header.textContent.trim() === headerText) {
                return header;
            }
        }
        return null;
    }
});
