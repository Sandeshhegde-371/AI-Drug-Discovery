import { BrowserRouter as Router, Routes, Route, Navigate } from 'react-router-dom';
import { OverviewPage } from './pages/OverviewPage';
import { UsesPage } from './pages/UsesPage';
import { SamplePage } from './pages/SamplePage';
import { HomePage } from './pages/HomePage';

export default function App() {
  return (
    <Router>
      <Routes>
        <Route path="/" element={<HomePage />} />
        <Route path="/overview" element={<OverviewPage />} />
        <Route path="/uses" element={<UsesPage />} />
        <Route path="/sample" element={<SamplePage />} />
        {/* Catch-all route - redirect to home */}
        <Route path="*" element={<Navigate to="/" replace />} />
      </Routes>
    </Router>
  );
}
