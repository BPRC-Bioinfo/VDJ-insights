document.addEventListener('DOMContentLoaded', function () {
    // Toggle additional content sections
    const showMoreButtons = document.querySelectorAll('.show-more-btn');
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
            const tdElements = table.querySelectorAll('td');

            thElements.forEach(th => {
                const header = th.innerText;
                const cellIndex = Array.from(th.parentElement.children).indexOf(th);
                if (selectedColumns.includes(header)) {
                    th.style.display = "";
                    table.querySelectorAll('tr').forEach(tr => {
                        if (tr.children[cellIndex]) {
                            tr.children[cellIndex].style.display = "";
                        }
                    });
                } else {
                    th.style.display = "none";
                    table.querySelectorAll('tr').forEach(tr => {
                        if (tr.children[cellIndex]) {
                            tr.children[cellIndex].style.display = "none";
                        }
                    });
                }
            });
        });
    }

    // Apply initial column visibility based on checked columns
    applyColumnVisibility();
});
