document.addEventListener('DOMContentLoaded', function () {
    const toggleBtn = document.getElementById('toggle-btn');
    const sidebar = document.getElementById('sidebar');
    const container = document.querySelector('.container');
    const dropdownToggles = document.querySelectorAll('.sidebar-content .dropdown-toggle');
    const logo = document.querySelector('.sidebar-header .logo');
    const ellipsisIcons = document.querySelectorAll('.sidebar-content .fa-ellipsis-h');

    if (toggleBtn && sidebar && container) {
        toggleBtn.addEventListener('click', function () {
            sidebar.classList.toggle('minimized');
            container.classList.toggle('shifted');

            if (sidebar.classList.contains('minimized')) {
                if (logo) logo.style.display = 'none'; // Hide the logo
                ellipsisIcons.forEach(icon => icon.style.display = 'none'); // Hide ellipsis icons
            } else {
                if (logo) logo.style.display = 'block'; // Show the logo
                ellipsisIcons.forEach(icon => icon.style.display = 'inline-block'); // Show ellipsis icons
            }

            dropdownToggles.forEach(function (toggle) {
                const parent = toggle.closest('.sidebar-dropdown');
                parent.classList.remove('active');
                toggle.setAttribute('aria-expanded', 'false');
            });
        });
    }

    dropdownToggles.forEach(function (toggle) {
        toggle.addEventListener('click', function (e) {
            e.preventDefault();
            const parent = this.closest('.sidebar-dropdown');
            const isActive = parent.classList.toggle('active');

            this.setAttribute('aria-expanded', isActive.toString());

            dropdownToggles.forEach(function (otherToggle) {
                const otherParent = otherToggle.closest('.sidebar-dropdown');
                if (otherParent !== parent) {
                    otherParent.classList.remove('active');
                    otherToggle.setAttribute('aria-expanded', 'false');
                }
            });
        });
    });

    window.addEventListener('click', function (event) {
        if (!event.target.closest('.sidebar-dropdown')) {
            document.querySelectorAll('.sidebar-content .sidebar-dropdown').forEach(function (dropdown) {
                dropdown.classList.remove('active');
                dropdown.querySelector('.dropdown-toggle').setAttribute('aria-expanded', 'false');
            });
        }
    });
});
