document.addEventListener('DOMContentLoaded', function() {
    const searchInput = document.getElementById('searchInput');
    const isMac = navigator.platform.toUpperCase().indexOf('MAC') >= 0;
    const shortcutIndicator = document.createElement('div');
    shortcutIndicator.classList.add('shortcut-indicator');
    shortcutIndicator.innerHTML = `<span class="shortcut-key">${isMac ? '⌘ K' : 'Ctrl K'}</span>`;
    
    // Insert the shortcut indicator into the DOM
    const searchContainer = document.querySelector('.search-container');
    searchContainer.appendChild(shortcutIndicator);

    // Search function for a given table by its ID and a specific column index
    function searchTable(tableId, columnIndex) {
        const table = document.getElementById(tableId);
        if (!table) {
            console.error(`Table with ID '${tableId}' not found.`);
            return;
        }

        const rows = table.getElementsByTagName('tr');

        // Loop through all table rows (except the header row)
        for (let i = 1; i < rows.length; i++) {
            let td = rows[i].getElementsByTagName('td')[columnIndex];
            if (td) {
                let textValue = td.textContent || td.innerText;
                if (textValue.toUpperCase().indexOf(searchInput.value.toUpperCase()) > -1) {
                    rows[i].style.display = ""; // Show row if it matches the search
                } else {
                    rows[i].style.display = "none"; // Hide row if no match
                }
            }
        }
    }

    // Function to search through all tables inside .library-card elements
    function searchTablesInCards() {
        const libraryCards = document.querySelectorAll('.library-card');

        libraryCards.forEach(card => {
            const rows = card.querySelectorAll('table.report-table tr');

            // Loop through all table rows (except the header row)
            let found = false; // To check if any row is found within this card
            for (let i = 1; i < rows.length; i++) {
                let td = rows[i].getElementsByTagName('td')[1]; // Searching in the second column (index 1)
                if (td) {
                    let textValue = td.textContent || td.innerText;
                    if (textValue.toUpperCase().indexOf(searchInput.value.toUpperCase()) > -1) {
                        rows[i].style.display = ""; // Show row if it matches the search
                        found = true;
                    } else {
                        rows[i].style.display = "none"; // Hide row if no match
                    }
                }
            }

            // Show or hide the entire card based on whether any row matches the search
            if (found) {
                card.style.display = "";
            } else {
                card.style.display = "none";
            }
        });
    }

    // Add event listener to the search input to filter both static tables and dynamic tables in cards
    searchInput.addEventListener('keyup', function() {
        if (document.getElementById('headerTable')) {
            searchTable('headerTable', 1);  // Search the second column in the header table
        }
        if (document.getElementById('annotationTable')) {
            searchTable('annotationTable', 1);  // Search the second column in the annotation table
        }
        searchTablesInCards();  // Search within dynamically generated tables in library cards
    });

    // Add event listener for Ctrl + K or Cmd + K to toggle focus on search input
    document.addEventListener('keydown', function(event) {
        if ((event.ctrlKey || event.metaKey) && event.key.toLowerCase() === 'k') {
            event.preventDefault(); // Prevent the default browser action
            
            if (document.activeElement === searchInput) {
                searchInput.blur(); // Remove focus if already focused
            } else {
                searchInput.focus(); // Focus on the search input if not already focused
            }
        }
    });
});
