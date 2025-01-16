document.addEventListener('DOMContentLoaded', function () {
    const toggleBtn = document.getElementById('toggle-btn');
    const toggleBtnHidden = document.getElementById('toggle-btn-hidden');
    const sidebar = document.getElementById('sidebar');
    const mainContent = document.getElementById('main-content');

    toggleBtn.addEventListener('click', function () {
        sidebar.classList.toggle('hidden');
        sidebar.classList.toggle('show');
        mainContent.classList.toggle('full-width');
        toggleBtn.classList.toggle('hidden');
        toggleBtnHidden.classList.toggle('hidden');
    });

    toggleBtnHidden.addEventListener('click', function () {
        sidebar.classList.toggle('hidden');
        sidebar.classList.toggle('show');
        mainContent.classList.toggle('full-width');
        toggleBtn.classList.toggle('hidden');
        toggleBtnHidden.classList.toggle('hidden');
    });
});
