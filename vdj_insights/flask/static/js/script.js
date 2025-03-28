const searchInput = document.getElementById('searchInput');
const amountTable = document.getElementById('amountTable');

const updateVisibleCount = () => {
    let visibleCount = 0;
    document.querySelectorAll('#DataTable tbody tr').forEach(row => {
        if (row.style.display !== 'none') {
            visibleCount++;
        }
    });
    amountTable.textContent = `Entries: ${visibleCount}`;
};

searchInput.addEventListener('keyup', function (event) {
    const searchTerm = event.target.value.toLowerCase();
    document.querySelectorAll('#DataTable tbody tr').forEach(row => {
        let rowText = '';
        row.querySelectorAll('td').forEach(td => {
            if (td.style.display !== 'none') {
                rowText += td.textContent.toLowerCase() + ' ';
            }
        });
        row.style.display = rowText.includes(searchTerm) ? '' : 'none';
    });
    updateVisibleCount();
});

document.addEventListener('DOMContentLoaded', updateVisibleCount);

document.getElementById('applyColumns').addEventListener('click', function () {
    let selectedColumns = [];
    document.querySelectorAll('.column-checkbox:checked').forEach(checkbox => {
        selectedColumns.push(checkbox.value);
    });

    document.querySelectorAll('.column-header').forEach(th => {
        const colName = th.getAttribute('data1-column');
        th.style.display = selectedColumns.includes(colName) ? "" : "none";
    });

    document.querySelectorAll('.column-data1').forEach(td => {
        const colName = td.getAttribute('data1-column');
        td.style.display = selectedColumns.includes(colName) ? "" : "none";
    });

    let modal = bootstrap.Modal.getInstance(document.getElementById('columnModal'));
    modal.hide();
});

document.getElementById('toggleSelectAll').addEventListener('click', function () {
    let checkboxes = document.querySelectorAll('.column-checkbox');
    let defaultColumns = [
        'Sample', 'Haplotype', 'Region', 'Segment',
        'Start coord', 'End coord', 'Reference', 'short name', 'Status'
    ];

    let allChecked = Array.from(checkboxes).every(checkbox => checkbox.checked);
    checkboxes.forEach(checkbox => {
        checkbox.checked = allChecked ? defaultColumns.includes(checkbox.value) : true;
    });

    this.textContent = allChecked ? 'Select All' : 'Reset';
});

document.addEventListener('DOMContentLoaded', function () {
    const table = document.getElementById('DataTable');
    const headers = table.querySelectorAll('thead th.sortable');

    headers.forEach((header, index) => {
        header.addEventListener('click', function () {
            const currentAsc = header.classList.contains('asc');
            sortTableByColumn(table, index, !currentAsc);

            headers.forEach(h => h.classList.remove('asc', 'desc'));
            if (!currentAsc) {
                header.classList.add('asc');
            } else {
                header.classList.add('desc');
            }
        });
    });

    function sortTableByColumn(table, columnIndex, ascending = true) {
        const tBody = table.querySelector('tbody');
        const rows = Array.from(tBody.querySelectorAll('tr'));

        rows.sort((rowA, rowB) => {
            const cellA = rowA.querySelectorAll('td')[columnIndex].textContent.trim();
            const cellB = rowB.querySelectorAll('td')[columnIndex].textContent.trim();

            const numA = parseFloat(cellA.replace(',', '.'));
            const numB = parseFloat(cellB.replace(',', '.'));
            const aIsNumeric = !isNaN(numA);
            const bIsNumeric = !isNaN(numB);

            if (aIsNumeric && bIsNumeric) {
                return ascending ? numA - numB : numB - numA;
            } else {
                if (cellA.toLowerCase() < cellB.toLowerCase()) r / eturn
                ascending ? -1 : 1;
                if (cellA.toLowerCase() > cellB.toLowerCase()) return ascending ? 1 : -1;
                return 0;
            }
        });

        rows.forEach(row => tBody.appendChild(row));
    }
});

document.addEventListener('DOMContentLoaded', function () {
    const searchInputs = document.querySelectorAll('input[id^="searchInput"]');

    searchInputs.forEach(searchInput => {
        const tableId = searchInput.id.replace('searchInput', 'DataTable');
        const amountButtonId = searchInput.id.replace('searchInput', 'amountTable');
        const table = document.getElementById(tableId);
        const amountButton = document.getElementById(amountButtonId);

        if (table && amountButton) {
            const tableRows = table.querySelectorAll('tbody tr');

            const updateVisibleCount = () => {
                let visibleCount = 0;

                tableRows.forEach(row => {
                    if (row.style.display !== 'none') {
                        visibleCount++;
                    }
                });

                amountButton.textContent = `Entries: ${visibleCount}`;
            };

            updateVisibleCount();

            searchInput.addEventListener('keyup', function (event) {
                const searchTerm = event.target.value.toLowerCase();

                tableRows.forEach(row => {
                    const rowText = row.textContent.toLowerCase();
                    if (rowText.includes(searchTerm)) {
                        row.style.display = '';
                    } else {
                        row.style.display = 'none';
                    }
                });

                updateVisibleCount();
            });
        } else {
            console.warn(`Table or button not found for search input "${searchInput.id}"`);
        }
    });
});