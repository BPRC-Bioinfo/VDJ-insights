document.addEventListener('DOMContentLoaded', function () {
    const toggleBtn = document.getElementById('toggle-btn');
    const dropdownToggles = document.querySelectorAll('.sidebar-content .dropdown-toggle');
    
    // Function to toggle the minimized state
    function toggleSidebar() {
        const isMinimized = document.documentElement.classList.toggle('sidebar-minimized');
        localStorage.setItem('sidebarMinimized', isMinimized);
    }

    if (toggleBtn) {
        toggleBtn.addEventListener('click', function () {
            toggleSidebar();

            // Reset dropdown menus
            dropdownToggles.forEach(function (toggle) {
                const parent = toggle.closest('.sidebar-dropdown');
                if (parent) {
                    parent.classList.remove('active');
                    toggle.setAttribute('aria-expanded', 'false');
                }
            });
        });
    }

    // Dropdown toggle functionality
    dropdownToggles.forEach(function (toggle) {
        toggle.addEventListener('click', function (e) {
            e.preventDefault();
            const parent = this.closest('.sidebar-dropdown');
            const isActive = parent.classList.toggle('active');

            this.setAttribute('aria-expanded', isActive.toString());

            // Close other dropdowns
            dropdownToggles.forEach(function (otherToggle) {
                const otherParent = otherToggle.closest('.sidebar-dropdown');
                if (otherParent !== parent) {
                    otherParent.classList.remove('active');
                    otherToggle.setAttribute('aria-expanded', 'false');
                }
            });
        });
    });

    // Close dropdowns when clicking outside
    window.addEventListener('click', function (event) {
        if (!event.target.closest('.sidebar-dropdown')) {
            document.querySelectorAll('.sidebar-content .sidebar-dropdown').forEach(function (dropdown) {
                dropdown.classList.remove('active');
                dropdown.querySelector('.dropdown-toggle').setAttribute('aria-expanded', 'false');
            });
        }
    });
});
